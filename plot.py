#!/usr/bin/env /usr/bin/python
import os,sys,re,gzip,csv,numpy
from math import fabs,fsum
from getopt import getopt
from random import random
from time import localtime,mktime,strptime,time
from calendar import timegm
from datetime import datetime
from subprocess import call,Popen,PIPE
from numpy import mean,std
from scipy import optimize

from ROOT import AddressOf
from ROOT import gROOT,gSystem,gStyle
from ROOT import TFile,TTree,TBranch,TGraph,TH1I,TH1F,TH2I,TH2F,TH3F,TProfile,THStack,TCanvas,TLegend,TGaxis,TIter

gSystem.Load('libAstro')
gSystem.Load('libMisc')
gSystem.Load('libEvent')
gSystem.Load('libDst')
from ROOT import THPKT1_DST_EVENT,THPKT1_DST_NOTICE,THPKT1_DST_TIME,THPKT1_DST_THRESHOLD,THPKT1_DST_VOLTS # must load libDst first

hs_regex = re.compile(r'(y\d{4}m\d{2}d\d{2}p\d{2})') # fully qualified hal stamp regular expression
HS_regex = re.compile(r'(y\d{4}m\d{2}d\d{2})p(\d{2})') # now capturing date stamp and part stamp value

class HalStamp:
  def __init__(self,hs=None):
    if hs != None:
      self.hs = hs
      date,part = re.match(HS_regex,hs).groups()
      self.dt = datetime(*strptime(date,'y%Ym%md%d')[:5])
      self.part = int(part)
    else:
      self.hs = None
      self.dt = None
      self.part = None

  def __repr__(self):
    if self.hs == None : return 'None'
    else : return self.hs

  def __str__(self):
    if self.hs == None : return 'None'
    else : return self.hs

  def __hash__(self):
    return int(timegm(self.dt.timetuple())) + self.part

  def __lt__(self,other): # x < y
    if self.dt == None or other.dt == None : return False
    if self.dt < other.dt : return True
    if self.dt == other.dt and self.part < other.part : return True
    else : return False

  def __le__(self,other): # x <= y
    if self.dt == None or other.dt == None : return False
    if self.dt < other.dt : return True
    if self.dt == other.dt and self.part <= other.part : return True
    else : return False

  def __eq__(self,other): # x == y
    if other == None:
      if self.dt == other : return True
      else : return False
    if self.dt == other.dt and self.part == other.part : return True
    else : return False

  def __ne__(self,other): # x != y
    if other == None:
      if self.dt != other : return True
      else : return False
    if self.dt != other.dt : return True
    if self.dt == other.dt and self.part != other.part : return True
    else : return False

  def __ge__(self,other): # x >= y
    if self.dt == None or other.dt == None : return False
    if self.dt > other.dt : return True
    if self.dt == other.dt and self.part >= other.part : return True
    else : return False

  def __gt__(self,other): # x > y
    if self.dt == None or other.dt == None : return False
    if self.dt > other.dt : return True
    if self.dt == other.dt and self.part > other.part : return True
    else : return False

def convert_time(packet,tevent=None):
  '''convert packet timestamp to (integral) epoch milliseconds'''
  if type(packet) == type(THPKT1_DST_NOTICE()):
    time_stamp = '%d %d %02d:%02d:%02d' % (packet.year,packet.day,packet.hour,packet.min,packet.sec)
    return int(1000*mktime(strptime(time_stamp,'%Y %j %H:%M:%S')) + packet.msec)
  elif type(packet) == type(THPKT1_DST_TIME()):
    hour = packet.sec/3600
    min = packet.sec%3600/60
    sec = packet.sec%3600%60
    time_stamp = '%d %d %02d:%02d:%02d' % (packet.year,packet.day,hour,min,sec)
    return int(1000*mktime(strptime(time_stamp,'%Y %j %H:%M:%S')) + packet.msec[tevent])
  else:
    return None

def find_nearest_tuple(listA,listB):
  '''listA and listB are lists of tuples of numbers and are sorted by the
  first number in their tuples.  For each tuple in listA, associate a tuple
  in listB whose first number is closest to the first number from the tuple
  of listA'''
  listC = []
  last_tupleB = listB[0]
  for tupleA in listA:
    last_diff = fabs(tupleA[0] - last_tupleB[0]) + 1
    new_listB = listB[listB.index(last_tupleB):]
    for tupleB in new_listB:
      if fabs(tupleA[0] - tupleB[0]) <= last_diff:
        last_tupleB = tupleB
        last_diff = fabs(tupleA[0] - last_tupleB[0])
      else:
        break
    listC.append(tupleA + last_tupleB)
  return listC

def compute_bins(array):
  if len(array) < 2 : return (3,array[0] - 0.1*array[0],array[0] + 0.1*array[0])
  elif len(array) == 0 : return(1,0,1)
  # interquartile range
  sarray = sorted(array)
  n = len(sarray)
  index = (n/2 + 1)/2;
  if n%4 == 0 or n%4 == 1:
    first = (sarray[index - 1] + sarray[index])/2.0
    third = (sarray[n - index - 1] + sarray[n - index])/2.0
  else:
    first = sarray[index - 1]
    third = sarray[n - index]
  # # bins, bin width
  min = sarray[0]
  max = sarray[-1] + 1
  binw = 2*(third - first)/n**(1/3.) # Freedman-Diaconis' choice
  if binw == 0 : nbins = 1
  else : nbins = int((max - min)/binw)
  return (nbins,min,max)

class Data:
  def __init__(self):
    self.files = {}
    self.data = {}
    self.get_LED()
    #self.get_noise_closed()
    self.get_hvcalib()

  # TODO: use method decorator
  def get_LED(self):
    find_command = 'find /home/findlay/data/middle_drum/ -type f -name "*led355*"'
    #find_command = 'find /home/findlay/data/middle_drum/ -type f -name "*y2008m11d07p02*led355*"'
    self.LED_data = self.get_data(find_command,kind='led355')
  def get_noise_closed(self):
    find_command = 'find /home/findlay/data/middle_drum/ -type f -name "*noise-closed*"'
    self.LED_data = self.get_data(find_command,kind='noise-closed')
  def get_hvcalib(self):
    find_command = 'find /home/findlay/data/middle_drum/ -type f -name "*hvcalib*"'
    self.LED_data = self.get_data(find_command,kind='hvcalib')

  def get_data(self,find_command,kind=''):
    find = Popen(find_command,stdout=PIPE,shell=True)
    file_names = re.split('\n',find.stdout.read())[:-1]
    find.stdout.close()
    self.files[kind] = {}
    for file_name in sorted(file_names):
      hal_stamp = HalStamp(re.search(hs_regex,file_name).group(0))
      self.files[kind][hal_stamp] = TFile(file_name)
    self.data[kind] = {}
    for hal_stamp in self.files[kind].keys():
      tree = self.files[kind][hal_stamp].Get('T')
      if tree == None : continue
      self.data[kind][hal_stamp] = {
          'time':tree.GetBranch('time'),
          'event':tree.GetBranch('event'),
          'snapshot':tree.GetBranch('snapshot'),
          'minute':tree.GetBranch('minute'),
          'threshold':tree.GetBranch('threshold'),
          'countrate':tree.GetBranch('countrate'),
          'volts':tree.GetBranch('volts'),
          'notice':tree.GetBranch('notice'),
          'remote':tree.GetBranch('remote'),
          'calib':tree.GetBranch('calib'),
          'boardid':tree.GetBranch('boardid'),
          'mstat':tree.GetBranch('mstat')}

  def print_packets(self,branch_type,kind='led355',hal_stamp=None):
    if hal_stamp == None:
      hal_stamps = sorted(self.data[kind].keys())
    else:
      hal_stamps = (hal_stamp,)
    for hal_stamp in hal_stamps:
      branch = self.data[kind][hal_stamp][branch_type]
      if branch_type == 'event'  : packet = THPKT1_DST_EVENT()
      if branch_type == 'notice' : packet = THPKT1_DST_NOTICE()
      if branch_type == 'time'   : packet = THPKT1_DST_TIME()
      if branch_type == 'volts'  : packet = THPKT1_DST_VOLTS()
      branch.SetAddress(AddressOf(packet))
      for entry in xrange(branch.GetEntries()):
        branch.GetEntry(entry)
        if branch_type == 'event'  : print hal_stamp,packet.pktHdr_crate,packet.event,packet.version,packet.minute,packet.msec,packet.ntubes
        if branch_type == 'notice' : print hal_stamp,packet.pktHdr_crate,packet.type,packet.year,packet.day,packet.hour,packet.min,packet.sec,packet.msec,'     ',packet.text
        if branch_type == 'time'   : print hal_stamp,packet.pktHdr_crate,packet.year,packet.day,packet.sec,packet.freq,packet.mark_error,packet.minute_offset,packet.error_flags,packet.events
        if branch_type == 'volts'  : print hal_stamp,packet.pktHdr_crate,packet.minute,packet.obVer,packet.hvChnls,packet.ob_p12v,packet.ob_p05v,packet.ob_n12v,packet.ob_n05v,packet.ob_tdcRef,packet.ob_temp,packet.ob_thRef,packet.ob_gnd,packet.garb_temp,packet.garb_p12v,packet.garb_n12v,packet.garb_p05v,packet.garb_s05v,packet.garb_lemo1,packet.garb_anlIn,packet.garb_clsVolts,packet.garb_clsTemp,packet.garb_mirX,packet.garb_mirY,packet.garb_clsX,packet.garb_clsY,packet.garb_ns,packet.garb_hvSup,packet.garb_hvChnl,packet.cluster,packet.hv
      print '%s  %5d' % (hal_stamp,branch.GetEntries())

class Plot:
  def __init__(self):
    gROOT.Reset()
    gROOT.SetStyle('Plain')
    gStyle.SetOptFit(11)
    gStyle.SetOptStat('nrme')

    self.data_object = Data()
    self.data = self.data_object.data
    if os.path.isfile('/home/findlay/data/flash_stats_cluster.gz'):
      self.flash_stats_cluster_file = gzip.open('/home/findlay/data/flash_stats_cluster.gz','r')
    if os.path.isfile('/home/findlay/data/flash_stats_tube.gz'):
      self.flash_stats_tube_file = gzip.open('/home/findlay/data/flash_stats_tube.gz','r')
    if os.path.isfile('/home/findlay/data/hv_calib.txt'):
      self.HV_calib_file = file('/home/findlay/data/hv_calib.txt','r')
    #hal_stamp := year month day part mirror
    self.flash_stats_cluster_regex = re.compile(r'^(\d+) (y\d{4}m\d{2}d\d{2}p\d{2})m(\d{2}) (\d+) (\d+) (\S+) (\S+) (\S+) (\S+) (\S+)$')
    self.flash_stats_tube_regex    = re.compile(r'^(\d+) (y\d{4}m\d{2}d\d{2}p\d{2})m(\d{2})t(\d{3}) (\d+) (\d+) (\S+) (\S+) (\S+) (\S+) (\S+)$')
    self.LED_temp_regex = re.compile(r'^TEMP A (\S+) B (\S+) C (\S+) D (\S+)$')
    self.PTH_regex = re.compile(r'^m(\d{1,2}): @\d+ (\S+) (\S+) (\S+)')
    self.canvas = TCanvas('LED','LED canvas',1024,791)

    self.t0 = time()
    self.t1 = self.diff_time(sorted(self.data['led355'].keys())[0])
    self.t2 = self.diff_time(sorted(self.data['led355'].keys())[-1])
    gStyle.SetTimeOffset(self.t0)

    self.time_bins = (6*(int((self.t2 - self.t1)/86400.) + 1),self.t1,self.t2 + 86400) # bin sizes are 6/day
    self.hent_bins = (384,0,384)
    self.mean_bins = (1024,0,2560)
    self.RMS_bins = (512,0,512)
    self.const_bins = (512,0,64)
    self.P_bins = (512,700,1000)
    self.T_bins = (1024,250,350)
    self.H_bins = (512,0,100)
    self.LED_AB_bins = (1024,318,318.5)
    self.LED_C_bins = (1024,317,318.5)
    self.volts_bins = (512,512,2048)
    self.peg_bins = (600,0,30)

    self.volts_ped = 1.97747874989

    self.mirrors = 14
    self.tubes = 256

  def get_Stan_HV_calib_data(self):
    self.HV_calib_stats_regex = re.compile('^(?P<mir>\d+) (?P<subcl>\d+) (?P<tube>\d+) (?P<m>\S+) (?P<b>\S+) (?P<m2>\S+) (?P<e>\S+) (?P<e2>\S+) (?P<n>\d+) (?P<n2>\d+) (?P<n0>\d+)$')
    self.HV_calib_stats = {}
    for line in self.HV_calib_file:
      m = re.match(self.HV_calib_stats_regex,line)
      if m:
        mirror = int(m.group('mir'))
        tube = 16*int(m.group('subcl')) + int(m.group('tube'))
        if mirror not in self.HV_calib_stats:
          self.HV_calib_stats[mirror] = {}
        if tube not in self.HV_calib_stats[mirror]:
          self.HV_calib_stats[mirror][tube] = {
              'm':float(m.group('m')),
              'b':float(m.group('b')),
              'm2':float(m.group('m2')),
              'e':float(m.group('e')),
              'e2':float(m.group('e2')),
              'n':int(m.group('n')),
              'n2':int(m.group('n2')),
              'n0':int(m.group('n0'))}

  def get_Doug_HV_pegs(self,file_name='/home/findlay/data/pegs/MD_Final_HV_Peggings.xls'):
    '''read Doug's spreadsheet'''
    antixls_command = ('antixls','--csv',file_name)
    antixls = Popen(antixls_command,stdout=PIPE)
    HV_pegs_csv = csv.reader(re.split('\n',antixls.stdout.read())[:-1])
    antixls.stdout.close()
    mirror = 0
    pegs = {}
    for line in HV_pegs_csv:
      if re.match(r'^Desired Value$',line[0]): # new header ==> new sheet (new mirror)
        mirror += 1
        pegs[mirror] = {}
      else:
        pegs[mirror][int(line[9])] = int(line[13]) # mirror.tube = peg
    return pegs

  def get_HV_pegs(self,file_name='/home/findlay/data/pegs/MD_2009-07-07_pegs.txt'):
    pegs = {}
    for mirror in xrange(1,self.mirrors + 1):
      pegs[mirror] = {}
    tube = -1
    for line in csv.reader(file(file_name),delimiter='\t'):
      if re.match(r'\d',line[0]):
        tube += 1
        for mirror in sorted(pegs.keys()):
          pegs[mirror][tube] = int(line[mirror + 2])
    return pegs

  def diff_time(self,hal_stamp):
    '''convert yYYmMMdDDpPP hal stamp into epoch time seconds and subtract t0'''
    #return mktime(strptime(re.match(r'(y\d{4}m\d{2}d\d{2})',hal_stamp.hs).group(0),'y%Ym%md%d')) - self.t0
    return mktime(hal_stamp.dt.timetuple()) - self.t0

  def get_run_start(self,hal_stamp,kind):
    notice_branch = self.data[kind][hal_stamp]['notice']
    notice = THPKT1_DST_NOTICE()
    notice_branch.SetAddress(AddressOf(notice))
    for entry in xrange(notice_branch.GetEntries()):
      notice_branch.GetEntry(entry)
      if notice.type == 8: # event time is measured as offset from first RUN START
        return convert_time(notice)

  def serialize_packets(self,hal_stamp,kind):
    start = self.get_run_start(hal_stamp,kind)

  def adj_stats_box(self,hist):
    stats = hist.GetListOfFunctions().FindObject('stats')
    stats.SetX1NDC(0.70) ; stats.SetX2NDC(0.90)
    stats.SetY1NDC(0.60) ; stats.SetY2NDC(0.90)
    stats.SetFillStyle(0)
    stats.SetTextFont(42)

  def adj_time_axis(self,hist):
    hist.GetXaxis().SetTimeDisplay(1)
    hist.GetXaxis().SetLabelSize(0.03)
    hist.GetXaxis().SetTimeFormat('%Y-%m-%d')

  def adj_hist_fonts(self,hist):
    hist.SetTitleFont(42)
    hist.SetLabelFont(42)

  def adj_axes(self,hist):
    hist.GetYaxis().SetTitleOffset(1.3)
    hist.GetXaxis().SetTitleFont(42)
    hist.GetXaxis().SetLabelFont(42)
    hist.GetYaxis().SetTitleFont(42)
    hist.GetYaxis().SetLabelFont(42)
    hist.GetZaxis().SetTitleFont(42)
    hist.GetZaxis().SetLabelFont(42)

  def adj_funcs(self,hist):
    next = TIter(hist.GetListOfFunctions())
    while True:
      try:
        i = next()
        i.SetLineWidth(2)
        i.SetNpx(1024)
      except:
        break

  def display_graph(self,graph,canvas=None,draw_option='a',stat_style='nrmeou',fit_style=11,pause=True):
    '''be sure to turn batch mode off'''
    if canvas == None:
      canvas = self.canvas
    canvas.cd()
    graph.Draw(draw_option)
    self.adj_axes(graph)
    self.adj_funcs(graph)
    gStyle.SetOptFit(fit_style)
    gStyle.SetOptStat(stat_style)
    if not fit_style == 0 or not stat_style == '':
      self.adj_stats_box(graph,canvas)
    canvas.Modified() ; canvas.Update()
    if pause == True:
      raw_input()

  def display_hist_stack(self,stack,canvas=None,legend=None,draw_option='',stat_style='nrmeou',fit_style=11,pause=True):
    '''be sure to turn batch mode off'''
    if canvas == None:
      canvas = self.canvas
    canvas.cd()
    stack.Draw(draw_option)
    hist = stack.GetHistogram()
    self.adj_axes(hist)
    self.adj_hist_fonts(hist)
    self.adj_funcs(hist)
    gStyle.SetOptFit(fit_style)
    gStyle.SetOptStat(stat_style)
#    if not fit_style == 0 or not stat_style == '':
#      self.adj_stats_box(hist,canvas)
    if legend != None:
      legend.Draw()
    canvas.SetTitle(hist.GetName())
    canvas.Modified() ; canvas.Update()
    if pause == True:
      raw_input()

  def display_histogram(self,hist,canvas=None,xaxis_time=False,draw_option='',stat_style='nrmeou',fit_style=11,pause=True):
    '''be sure to turn batch mode off'''
    if canvas == None:
      canvas = self.canvas
    canvas.cd()
    hist.Draw(draw_option)
    self.adj_axes(hist)
    self.adj_hist_fonts(hist)
    self.adj_funcs(hist)
    gStyle.SetOptFit(fit_style)
    gStyle.SetOptStat(stat_style)
    if not fit_style == 0 or not stat_style == '':
      self.adj_stats_box(hist)
    if xaxis_time:
      self.adj_time_axis(hist)
    canvas.SetTitle(hist.GetName())
    canvas.Modified() ; canvas.Update()
    if pause == True:
      raw_input()

  def write_kind_plots(self,hist_dict,plot_file,xaxis_time=False,draw_option='',stat_style='nrmeou',fit_style=11,convert_file=True):
    '''accepts a dict of dicts of hists.  For each dict in the main dict, create a multipage plot file'''
    hists = {}
    for mirror_hists in hist_dict.values():
      for hist_kind in mirror_hists.keys():
        if not hist_kind in hists:
          hists[hist_kind] = [mirror_hists[hist_kind]]
        else:
          hists[hist_kind].append(mirror_hists[hist_kind])
    for kind in hists.keys():
      self.write_plots(plot_file % kind,hists[kind],xaxis_time,draw_option)
    if convert_file == True:
      self.convert_plot_file(plot_file)

  def write_plots(self,plot_file,hist_list,xaxis_time=False,draw_option='',stat_style='nrmeou',fit_style=11,convert_file=True):
    self.canvas.cd()
    self.canvas.Print('%s[' % plot_file)
    for hist in hist_list:
      self.write_plot(plot_file,hist,xaxis_time,draw_option,stat_style,fit_style,None,False)
    self.canvas.Print('%s]' % plot_file)
    if convert_file == True:
      self.convert_plot_file(plot_file)

  def write_plot(self,plot_file,hist,canvas=None,xaxis_time=False,draw_option='',stat_style='nrmeou',fit_style=11,legend=None,convert_file=True):
    if canvas == None:
      canvas = self.canvas
    hist.Draw(draw_option)
    self.canvas.Modified() ; self.canvas.Update()
    if legend == None: # THStacks are strange
      self.adj_axes(hist)
      self.adj_hist_fonts(hist)
      self.adj_funcs(hist)
      gStyle.SetOptFit(fit_style)
      gStyle.SetOptStat(stat_style)
      if not fit_style == 0 or not stat_style == '':
        self.adj_stats_box(hist,canvas)
    else:
      legend.Draw()
    if xaxis_time:
      self.adj_time_axis(hist)
    self.canvas.Modified() ; self.canvas.Update()
    if not os.path.isdir(os.path.split(plot_file)[0]):
      os.makedirs(os.path.split(plot_file)[0])
    self.canvas.Print(plot_file)

  def convert_plot_file(self,plot_file):
    base,ps_file = os.path.split(plot_file)
    pdf_file = os.path.splitext(ps_file)[0] + '.pdf'
    os.chdir(base)
    call(('ps2pdf',ps_file,pdf_file))
#    if os.path.isfile('/tmp/gs_*'): os.remove('/tmp/gs_*')
    print '%s/%s' % (base,pdf_file)
