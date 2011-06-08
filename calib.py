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

sys.argv.append('-b') # run in batch mode
if '-d' in sys.argv : del sys.argv[sys.argv.index('-b')] # run in display mode

from ROOT import AddressOf
from ROOT import gROOT,gSystem,gStyle
from ROOT import TFile,TTree,TBranch,TGraph,TH1I,TH1F,TH2I,TH2F,TH3F,TProfile,THStack,TCanvas,TLegend,TGaxis,TIter

gSystem.Load('libAstro')
gSystem.Load('libMisc')
gSystem.Load('libEvent')
gSystem.Load('libDst')
from ROOT import THPKT1_DST_EVENT,THPKT1_DST_NOTICE,THPKT1_DST_TIME,THPKT1_DST_THRESHOLD,THPKT1_DST_VOLTS # must load libDst first

from plot import HalStamp,Plot
from plot import convert_time,find_nearest_tuple,compute_bins

class HVCalib(Plot):
  '''investigate HV calibrations, and problems with the HV system'''

  def get_HV_calib_steps(self,hal_stamp):
    '''return voltage steps used as the calibration voltage supply'''
    step_voltages = {}
    if hal_stamp == 'y2009m02d14p01' or hal_stamp == 'y2009m03d01p01':
      '''iterate through the notice packets to find the manually entered supply voltages, and structure them by their associated mirror'''
      notice_branch = self.data['hvcalib'][hal_stamp]['notice']
      notice = THPKT1_DST_NOTICE()
      notice_branch.SetAddress(AddressOf(notice))
      mirror = 1
      for entry in xrange(notice_branch.GetEntries()):
        notice_branch.GetEntry(entry)
        try:
          voltage_step = float(re.search(r'^([0-9.]+)$',notice.text).group(1))
        except:
          continue
        if len(step_voltages) == 0:
          step_voltages[mirror] = [voltage_step]
        elif voltage_step < step_voltages[mirror][-1]:
          mirror += 1
          step_voltages[mirror] = [voltage_step]
        else:
          step_voltages[mirror].append(voltage_step)
    elif hal_stamp == 'y2009m08d10p01': pass
    elif hal_stamp == 'y2009m08d10p02': pass
    elif hal_stamp == 'y2009m08d10p03': pass
    elif hal_stamp == 'y2009m08d10p04': pass
    elif hal_stamp == 'y2009m08d10p05': pass
    elif hal_stamp == 'y2009m08d10p06': pass
    elif hal_stamp == 'y2009m08d10p07': pass
    elif hal_stamp == 'y2009m08d10p08': pass
    elif hal_stamp == 'y2009m08d10p09': pass
    elif hal_stamp == 'y2009m08d10p10': pass
    elif hal_stamp == 'y2009m08d10p11': pass
    elif hal_stamp == 'y2009m08d10p12': pass
    elif hal_stamp == 'y2009m08d21p01':
      for mirror in xrange(1,self.mirrors + 1):
        step_voltages[mirror] = (401,)
    elif hal_stamp == 'y2009m08d21p02':
      for mirror in xrange(1,self.mirrors + 1):
        step_voltages[mirror] = (597,)
    elif hal_stamp == 'y2009m08d21p03':
      for mirror in xrange(1,self.mirrors + 1):
        step_voltages[mirror] = (596,)
    elif hal_stamp == 'y2009m08d21p04':
      for mirror in xrange(1,self.mirrors + 1):
        step_voltages[mirror] = (796,)
    elif hal_stamp == 'y2009m08d21p05':
      for mirror in xrange(1,self.mirrors + 1):
        step_voltages[mirror] = (1000,)
    elif hal_stamp == 'y2009m08d21p06':
      for mirror in xrange(1,self.mirrors + 1):
        step_voltages[mirror] = (1197,)
    elif hal_stamp == 'y2009m08d21p07':
      for mirror in xrange(1,self.mirrors + 1):
        step_voltages[mirror] = (1398,)
    return step_voltages

  def get_HV_calib(self):
    HV_calib = {}
    for hal_stamp in self.data['hvcalib'].keys():
      if re.match(r'y2009m08d10',hal_stamp) : continue
      HV_calib[hal_stamp] = {}
      step_voltages = self.get_HV_calib_steps(hal_stamp)
      volts_branch = self.data['hvcalib'][hal_stamp]['volts']
      volts = THPKT1_DST_VOLTS()
      volts_branch.SetAddress(AddressOf(volts))
      for mirror in sorted(step_voltages.keys()):
        if not mirror == 6 : continue
        HV_calib[hal_stamp][mirror] = {}
        for entry in xrange(volts_branch.GetEntries()):
          volts_branch.GetEntry(entry)
          if not mirror == volts.pktHdr_crate : continue
          for tube in xrange(volts.hvChnls):
            mean = float(sum(volts.hv))/len(volts.hv)
            step_voltage = min(step_voltages[mirror],key=lambda step : fabs(step - mean))
            if not tube in HV_calib[hal_stamp][mirror]:
              HV_calib[hal_stamp][mirror][tube] = {'measured volts':[],'supplied volts':[]}
            HV_calib[hal_stamp][mirror][tube]['measured volts'].append(volts.hv[tube] - self.volts_ped)
            HV_calib[hal_stamp][mirror][tube]['supplied volts'].append(step_voltage)
    return HV_calib

  def get_HV_calib_pars(self,fit='2par'):
    HV_calib = self.get_HV_calib()
    HV_calib_pars = {}
    for hal_stamp in HV_calib:
      HV_calib_pars[hal_stamp] = {}
      for mirror in HV_calib[hal_stamp]:
        HV_calib_pars[hal_stamp][mirror] = {}
        for tube in HV_calib[hal_stamp][mirror]:
          x = HV_calib[hal_stamp][mirror][tube]['measured volts']
          y = HV_calib[hal_stamp][mirror][tube]['supplied volts']
          if fit == '1par':
            fitfunc = lambda p,x : p[0]*x
            errfunc = lambda p,x,y : fitfunc(p,x) - y
            p0 = [1.] # parameter initial values
            p1,covariance,details,mesg,success = optimize.leastsq(errfunc,p0,args=(x,y),full_output=True)
            HV_calib_pars[hal_stamp][mirror][tube] = {'m':p1[0],'m_err':covariance[0][0]}
          elif fit == '2par':
            fitfunc = lambda p,x : p[1]*x + p[0]
            errfunc = lambda p,x,y : fitfunc(p,x) - y
            p0 = [0.,1]
            p1,covariance,details,mesg,success = optimize.leastsq(errfunc,p0,args=(x,y),full_output=True)
            HV_calib_pars[hal_stamp][mirror][tube] = {'m':p1[1],'m_err':covariance[1][1],'b':p1[0],'b_err':covariance[0][0]}
    return HV_calib_pars

  def plot_HV_calib_pars(self):
    HV_calib_pars = self.get_HV_calib_pars()

    m = []
    m_hist_stack = THStack('fitted slopes','HV calibration fitted slopes')
    ml = TLegend(0.80,0.50,0.90,0.90)
    for mirror in HV_calib_pars:
      m_hist = TH1F('m%02d slope' % mirror,';m%02d fit slope' % mirror,128,0.8,1.2)
      m_hist.SetFillColor(mirror)
      m_hist_stack.Add(m_hist)
      ml.AddEntry(m_hist,'m%02d' % mirror,'f')
      for tube in HV_calib_pars[mirror]:
        m.append(HV_calib_pars[mirror][tube]['m'])
        m_hist.Fill(HV_calib_pars[mirror][tube]['m'])
    print compute_bins(m)
    self.display_hist_stack(m_hist_stack,legend=ml)

#    m_err = []
#    m_err_hist_stack = THStack('fitted slopes','HV calibration fitted slopes errors')
#    m_errl = TLegend(0.80,0.50,0.90,0.90)
#    for mirror in HV_calib_pars:
#      m_err_hist = TH1F('m%02d slope error' % mirror,';m%02d fit slope error' % mirror,128,0,5)
#      m_err_hist.SetFillColor(mirror)
#      m_err_hist_stack.Add(m_err_hist)
#      m_errl.AddEntry(m_err_hist,'m%02d' % mirror,'f')
#      for tube in HV_calib_pars[mirror]:
#        m_err.append(HV_calib_pars[mirror][tube]['m_err'])
#        m_err_hist.Fill(HV_calib_pars[mirror][tube]['m_err'])
#    print compute_bins(m_err)
#    self.display_hist_stack(m_err_hist_stack,legend=m_errl)

    b = []
    b_hist_stack = THStack('fitted intercepts','HV calibration fitted intercepts')
    bl = TLegend(0.80,0.50,0.90,0.90)
    for mirror in HV_calib_pars:
      b_hist = TH1F('m%02d intercept' % mirror,';m%02d fit intercept' % mirror,256,-13.082888426084153,613.57714578208049)#-20,50)
      b_hist.SetFillColor(mirror)
      b_hist_stack.Add(b_hist)
      bl.AddEntry(b_hist,'m%02d' % mirror,'f')
      for tube in HV_calib_pars[mirror]:
        b.append(HV_calib_pars[mirror][tube]['b'])
        b_hist.Fill(HV_calib_pars[mirror][tube]['b'])
    print compute_bins(b)
    self.display_hist_stack(b_hist_stack,legend=bl)

  def plot_HV_calib_tubes(self):
    HV_calib_pars = self.get_HV_calib_pars()
    for mirror in HV_calib_pars:
      raw_volts = {}
      for hal_stamp in sorted(self.data['led355'].keys()):
        volts_branch = self.data['led355'][hal_stamp]['volts']
        volts = THPKT1_DST_VOLTS()
        volts_branch.SetAddress(AddressOf(volts))
        for entry in xrange(volts_branch.GetEntries()):
          volts_branch.GetEntry(entry)
          if not mirror == volts.pktHdr_crate : continue
          for tube in xrange(volts.hvChnls):
            if not tube in raw_volts:
              raw_volts[tube] = []
            raw_volts[tube].append(volts.hv[tube])
      for tube in raw_volts:
        mean = float(sum(raw_volts[tube]))/len(raw_volts[tube])
        hist = TH1F('m%02dt%03d' % (mirror,tube),';m%02dt%03d volts [V]' % (mirror,tube),1600,700,1500)
        vals = []
        for val in raw_volts[tube]:
          hist.Fill(val)
          vals.append(val)
        print mirror,tube,vals,float(sum(raw_volts[tube]))/len(raw_volts[tube])
        self.display_histogram(hist)

  def plot_calibrated_HV(self):
    HV_pegs = self.get_HV_pegs()
    HV_calib_pars = self.get_HV_calib_pars()
    plot_file = '/home/findlay/data/plots/HV/calibrated_HV.ps'
    self.canvas.Print('%s[' % plot_file)
    for hal_stamp in HV_calib_pars:
      for mirror in HV_calib_pars[hal_stamp]:
        HV = {}
        for hal_stamp in sorted(self.data['led355'].keys()):
          volts_branch = self.data['led355'][hal_stamp]['volts']
          volts = THPKT1_DST_VOLTS()
          volts_branch.SetAddress(AddressOf(volts))
          for entry in xrange(volts_branch.GetEntries()):
            volts_branch.GetEntry(entry)
            if not mirror == volts.pktHdr_crate : continue
            for tube in xrange(volts.hvChnls):
              if not tube in HV:
                HV[tube] = []
              #voltage = volts.hv[tube] + random() - self.volts_ped
              voltage = HV_calib_pars[hal_stamp][mirror][tube]['m']*(volts.hv[tube] + random() - self.volts_ped)
              HV[tube].append(voltage)
        name = 'm%02d volts' % mirror
        title = ';m%02d pegs;m%02d volts [V]' % (mirror,mirror)
        hvcalib_hist = TH2F(name,title,*(self.peg_bins + self.volts_bins))
        hvcalib_hist.SetMarkerStyle(2)
        for tube in HV:
          mean = float(sum(HV[tube]))/len(HV[tube])
          hvcalib_hist.Fill(HV_pegs[mirror][tube],mean)
        hvcalib_hist.Fit('pol1','Q')
        self.display_histogram(hvcalib_hist,stat_style='nrme')
        self.write_plot(plot_file,hvcalib_hist,stat_style='n')
    self.canvas.Print('%s]' % plot_file)
    self.convert_plot_file(plot_file)

  def plot_HV_supply_vs_time(self):
    HV_supply_hists = {}
    for hal_stamp in sorted(self.data['led355'].keys()):
      start = self.get_run_start(hal_stamp,'led355')
      volts_branch = self.data['led355'][hal_stamp]['volts']
      volts = THPKT1_DST_VOLTS()
      volts_branch.SetAddress(AddressOf(volts))
      for entry in xrange(volts_branch.GetEntries()):
        volts_branch.GetEntry(entry)
        time_diff = (start + 60*1000*volts.minute)/1000. - self.t0
        mirror = volts.pktHdr_crate
        if not mirror in HV_supply_hists:
          HV_supply_hists[mirror] = TH2F('m%02d supply volts' % mirror,';time (LED runs);m%02d garb HV supply (lemo1) [V]' % mirror,*(self.time_bins + (2048,1200,1500)))
        HV_supply_hists[mirror].Fill(time_diff,volts.garb_lemo1 + random())
    self.write_plots('/home/findlay/data/plots/HV/HV_supply_volts.ps',HV_supply_hists.values(),xaxis_time=True,convert_file=True)

  def plot_HV_vs_time(self):
    plot_file = '/home/findlay/data/plots/HV/HV_vs_time.ps'
    self.canvas.Print('%s[' % plot_file)
    for mirror in xrange(1,self.mirrors + 1):
      if not mirror == 1 : continue
      for tube in xrange(self.tubes):
        h = TH2F('m%02dt%03d supply volts' % (mirror,tube),';time (LED runs);m%02dt%03d HV [V]' % (mirror,tube),*(self.time_bins + (2048,512,2048)))
        for hal_stamp in sorted(self.data['led355'].keys()):
          start = self.get_run_start(hal_stamp,'led355')
          volts_branch = self.data['led355'][hal_stamp]['volts']
          volts = THPKT1_DST_VOLTS()
          volts_branch.SetAddress(AddressOf(volts))
          for entry in xrange(volts_branch.GetEntries()):
            volts_branch.GetEntry(entry)
            if not volts.pktHdr_crate == mirror : continue
            time_diff = (start + 60*1000*volts.minute)/1000. - self.t0
            h.Fill(time_diff,volts.hv[tube] + random())
        self.display_histogram(h,xaxis_time=True)
        self.write_plot(plot_file,h,xaxis_time=True)
    self.canvas.Print('%s]' % plot_file)
    self.convert_plot_file(plot_file)

  def plot_HV_supply_vs_time(self):
    HV_supply_hists = {}
    for mirror in xrange(1,self.mirrors + 1):
      HV_supply_hists[mirror] = TH2F('m%02d HV supply volts' % mirror,';time (LED runs);m%02d HV [V]' % mirror,*(self.time_bins + (2048,512,2048)))
      for hal_stamp in sorted(self.data['led355'].keys()):
        start = self.get_run_start(hal_stamp,'led355')
        volts_branch = self.data['led355'][hal_stamp]['volts']
        volts = THPKT1_DST_VOLTS()
        volts_branch.SetAddress(AddressOf(volts))
        for entry in xrange(volts_branch.GetEntries()):
          volts_branch.GetEntry(entry)
          if not volts.pktHdr_crate == mirror : continue
          time_diff = (start + 60*1000*volts.minute)/1000. - self.t0
          HV_supply_hists[mirror].Fill(time_diff,volts.garb_lemo1 + random())
    self.write_plots('/home/findlay/data/plots/HV/HV_supply_vs_time.ps',HV_supply_hists.values(),xaxis_time=True)

  def plot_HV_sub_vs_time(self,sub_type='subcluster'):
    HV_sub_hists = {}
    for sub in xrange(1,17):
      HV_sub_hists[sub] = TH2F('m06 subcl %d HV' % sub,';time (LED runs);%s %d HV [V]' % (sub_type,sub),*(self.time_bins + (2048,512,2048)))
      for hal_stamp in sorted(self.data['led355'].keys()):
        start = self.get_run_start(hal_stamp,'led355')
        volts_branch = self.data['led355'][hal_stamp]['volts']
        volts = THPKT1_DST_VOLTS()
        volts_branch.SetAddress(AddressOf(volts))
        for entry in xrange(volts_branch.GetEntries()):
          volts_branch.GetEntry(entry)
          if not volts.pktHdr_crate == 6 : continue
          time_diff = (start + 60*1000*volts.minute)/1000. - self.t0
          for tube in xrange(self.tubes):
            if sub_type == 'subcluster':
              if tube/16 + 1 == sub:
                HV_sub_hists[sub].Fill(time_diff,volts.hv[tube] + random())
            if sub_type == 'subtube':
              if tube%16 + 1 == sub:
                HV_sub_hists[sub].Fill(time_diff,volts.hv[tube] + random())
    self.write_plots('/home/findlay/data/plots/HV/HV_%s_vs_time.ps' % sub_type,HV_sub_hists.values(),xaxis_time=True)

  def plot_m06_20090821_HV_calib(self):
    hal_stamps = sorted(self.data['hvcalib'].keys())
    for hal_stamp in sorted(self.data['hvcalib'].keys()):
      if not re.search(r'y2009m08d21',hal_stamp) : continue
      mirror = 6
      step_voltages = self.get_HV_calib_steps(hal_stamp)
      volts_branch = self.data['hvcalib'][hal_stamp]['volts']
      volts = THPKT1_DST_VOLTS()
      volts_branch.SetAddress(AddressOf(volts))
      for entry in xrange(volts_branch.GetEntries()):
        volts_branch.GetEntry(entry)
        if not volts.pktHdr_crate == mirror : continue
        for tube in xrange(volts.hvChnls):
          if not tube in raw_volts:
            raw_volts[tube] = []
          raw_volts[tube].append(volts.hv[tube])

class LEDFlashes(Plot):
  '''plot each LED flash from the led355 data parts, find the mean,stddev,etc. of these plots'''

  def plot_QDCB_vs_time(self):
    for mirror in xrange(1,self.mirrors + 1):
      for tube in xrange(self.tubes):
        for hal_stamp in sorted(self.data['led355'].keys()):
          event_branch = self.data['led355'][hal_stamp]['event']
          event = THPKT1_DST_EVENT()
          event_branch.SetAddress(AddressOf(event))
          h = TProfile('%sm%02dt%03d QDCB' % (hal_stamp,mirror,tube),';event;%sm%02dt%03d QDCB' % (hal_stamp,mirror,tube),*(512,0,512))
          for entry in xrange(event_branch.GetEntries()):
            event_branch.GetEntry(entry)
            if mirror != event.pktHdr_crate : continue
            if tube in event.tube_num:
              h.Fill(event.event,event.qdcB[tube])
              h.Fit('pol1','Q')
#              h.Fit('pol0','QA')
          self.display_histogram(h)

  def observe_low_flashes(self):
    flash_times = {}
    for hal_stamp in sorted(self.data['led355'].keys()):
      flash_times[hal_stamp] = {}
      start = self.get_run_start(hal_stamp,'led355')
      time_branch = self.data['led355'][hal_stamp]['time']
      time = THPKT1_DST_TIME()
      time_branch.SetAddress(AddressOf(time))
      event_branch = self.data['led355'][hal_stamp]['event']
      event = THPKT1_DST_EVENT()
      event_branch.SetAddress(AddressOf(event))
      for entry in xrange(time_branch.GetEntries()):
        time_branch.GetEntry(entry)
        # get corresponding event for this tevent
        last_event_entry = 0
        for tevent in xrange(time.events):
          for event_entry in xrange(last_event_entry,event_branch.GetEntries()):
            event_branch.GetEntry(event_entry)
#            print start + 60*1000*event.minute + event.msec - convert_time(time,tevent),entry,event_entry,event.event ; raw_input()
            if start + 60*1000*event.minute + event.msec == convert_time(time,tevent):
              last_event_entry = event_entry
              break
          mirror = time.mirror[tevent]
          if not mirror in flash_times[hal_stamp]:
            flash_times[hal_stamp][mirror] = {'mean':0,'entries':[]}
          flash_times[hal_stamp][mirror]['entries'].append(time.nsec[tevent]%5e7)
        flash_times[hal_stamp][mirror]['mean'] = mean([t%5e7 for t in time.nsec[tevent]])
      print '%s  %5d' % (hal_stamp,time_branch.GetEntries())
    flash_hist = TH1F('combined detector nsecs',';nanoseconds%5e7 (normalized to 50 ms) [s];events',10000,0,1e8)
    for part in flash_times.values():
      for mirror in part.keys():
        for entry in part[mirror]['entries']:
          flash_hist.Fill(entry + (5e7 - part[mirror]['mean'])) # normalize all data means to 5e7
    self.canvas.cd()
    self.canvas.SetLogy()
    self.write_plot('/home/findlay/data/plots/LED/LED_resonance.ps',flash_hist,xaxis_time=False)

  def observe_flashes(self):
    flash_hists = {}
    for hal_stamp in self.data['led355'].keys():
      start = self.get_run_start(hal_stamp,'led355')
      event_branch = self.data['led355'][hal_stamp]['event']
      event = THPKT1_DST_EVENT()
      event_branch.SetAddress(AddressOf(event))
      for entry in xrange(event_branch.GetEntries()):
        event_branch.GetEntry(entry)
        h = TH1I('%sm%de%f' % (hal_stamp,event.pktHdr_crate,entry),';QDCB;tubes',*(compute_bins(event.qdcB)))
        for qdcB in event.qdcB : h.Fill(qdcB + random())
        h.Fit('gaus','LL Q')
        if h.GetMean() > 800:
          if not event.pktHdr_crate in flash_hists:
            flash_hists[event.pktHdr_crate] = h
        if len(flash_hists) == 14 : break
      print '%s  %5d' % (hal_stamp,event_branch.GetEntries())
    self.write_plots('/home/findlay/data/plots/LED/LED_flashes.ps',flash_hists.values(),xaxis_time=False)

  def plot_thresholds(self):
    threshold_hists = {}
    for hal_stamp in self.data['led355'].keys():
      start = self.get_run_start(hal_stamp,'led355')
      threshold_branch = self.data['led355'][hal_stamp]['threshold']
      threshold = THPKT1_DST_THRESHOLD()
      threshold_branch.SetAddress(AddressOf(threshold))
      for entry in xrange(threshold_branch.GetEntries()):
        threshold_branch.GetEntry(entry)
        mirror = threshold.pktHdr_crate
        h = TH1I('%sm%de%d' % (hal_stamp,mirror,entry),';threshold;tubes',*(compute_bins(threshold.thB)))
        for thB in threshold.thB : h.Fill(thB + random())
        h.Fit('gaus','LL Q')
        self.adj_stats_box(h,self.canvas)
        h.Draw()
        self.canvas.Modified() ; self.canvas.Update()
        raw_input()
        if not mirror in threshold_hists:
          threshold_hists[mirror] = TH2F('m%02d thresholds' % mirror,';time [200[8|9]-mm-dd];m%02d mean thB' % mirror,*(self.time_bins + (1000,0,1000)))
        time_diff = (start + 60*1000*threshold.min)/1000. - self.t0
        threshold_hists[mirror].Fill(time_diff,h.GetMean())
    self.write_plots('/home/findlay/data/plots/LED/LED_thresholds.ps',threshold_hists.values(),xaxis_time=True)

  def compute_flash_stats_cluster(self,flash_stats_cluster_file,hal_stamps,i):
    for hal_stamp in hal_stamps:
      start = self.get_run_start(hal_stamp,'led355')
      event_branch = self.data['led355'][hal_stamp]['event']
      event = THPKT1_DST_EVENT()
      event_branch.SetAddress(AddressOf(event))
      for entry in xrange(event_branch.GetEntries()):
        event_branch.GetEntry(entry) # root copies each object into 'event'
        h = TH1I('%sm%02d %d' % (hal_stamp,event.pktHdr_crate,i),';QDCB;tubes',*(compute_bins(event.qdcB)))
        # QDC truncates value to integer.  This spreads out the data randomly
        # throughout the interval like it originally was
        for qdcB in event.qdcB : h.Fill(qdcB + random())
        #h.Fit('gaus','E I LL M Q')
        h.Fit('gaus','LL Q')
        gaus = h.GetListOfFunctions().FindObject('gaus')
        hents,hmean,hRMS = h.GetEntries(),h.GetMean(),h.GetRMS()
        const,mean,sigma = gaus.GetParameter(0),gaus.GetParameter(1),gaus.GetParameter(2)
        # event time is measured as offset from first RUN START in milliseconds
        t = start + 60*1000*event.minute + event.msec
        flash_stats_cluster_file.write('%d %sm%02d %d %d %f %f %f %f %f\n' % (i,hal_stamp,event.pktHdr_crate,t,hents,hmean,hRMS,const,mean,sigma))
        i += 1
      print '%s  %5d' % (hal_stamp,event_branch.GetEntries())

  def compute_flash_stats_tube(self,flash_stats_file,hal_stamps,i):
    for hal_stamp in hal_stamps:
      events = {}
      start = self.get_run_start(hal_stamp,'led355')
      event_branch = self.data['led355'][hal_stamp]['event']
      event = THPKT1_DST_EVENT()
      event_branch.SetAddress(AddressOf(event))
      for entry in xrange(event_branch.GetEntries()):
        event_branch.GetEntry(entry)
        mirror = event.pktHdr_crate
        for i in xrange(event.ntubes):
          tube = event.tube_num[i]
          if not mirror in events:
            events[mirror] = {}
          if not tube in events[mirror]:
            events[mirror][tube] = []
          events[mirror][tube].append((event.qdcB[i],start + 60*1000*event.minute + event.msec))
      for mirror in events:
        for tube in events[mirror]:
          h = TH1I('%d' % i,';QDCB;events',*(compute_bins(events[mirror][tube][0])))
          times = []
          for QDCB,time_stamp in events[mirror][tube]:
            h.Fill(QDCB + random())
            times.append(time_stamp)
          h.Fit('gaus','LL Q')
          gaus = h.GetListOfFunctions().FindObject('gaus')
          hents,hmean,hRMS = h.GetEntries(),h.GetMean(),h.GetRMS()
          const,mean,sigma = gaus.GetParameter(0),gaus.GetParameter(1),gaus.GetParameter(2)
          t = sum(times)/len(times) # average timestamp
          flash_stats_file.write('%d %sm%02dt%03d %d %d %f %f %f %f %f\n' % (i,hal_stamp,mirror,tube,t,hents,hmean,hRMS,const,mean,sigma))
          i += 1
      print '%s  %5d' % (hal_stamp,event_branch.GetEntries())

  def make_flash_stats(self,cluster_file,tube_file):
    for file_name in (cluster_file,tube_file):
      if not os.path.isfile(file_name): # new file
        if not os.path.isdir(os.path.split(file_name)[0]):
          os.makedirs(os.path.split(file_name)[0])
        i = 0
        hal_stamps = sorted(self.data['led355'].keys())
      else: # append to existing file
        zcat = Popen(('/bin/zcat',file_name),stdout=PIPE)
        tail = Popen(('/usr/bin/tail','-n','1'),stdin=zcat.stdout,stdout=PIPE)
        #line = re.split('\n',tail.stdout.read())[0]
        line = tail.stdout.read()
        zcat.stdout.close() ; tail.stdout.close()
        if file_name == cluster_file : groups = re.match(self.flash_stats_cluster_regex,line).groups()
        elif file_name == tube_file : groups = re.match(self.flash_stats_tube_regex,line).groups()
        i,hal_stamp = int(groups[0]) + 1,HalStamp(groups[1])
        sorted_hal_stamps = sorted(self.data['led355'].keys())
        hal_stamps = sorted_hal_stamps[sorted_hal_stamps.index(hal_stamp) + 1:]
      flash_stats_file = gzip.open(file_name,'a+',9)
      if file_name == cluster_file : self.compute_flash_stats_cluster(flash_stats_file,hal_stamps,i)
      elif file_name == tube_file : self.compute_flash_stats_tube(flash_stats_file,hal_stamps,i)
      flash_stats_file.close()

class LEDEnv(Plot):
  '''plot LED flash mean,stddev,etc. versus various environmental parameters'''

  def get_flash_stats_cluster(self):
    stat_tuples = {}
    for line in self.flash_stats_cluster_file:
      match = re.match(self.flash_stats_cluster_regex,line).groups()
      i,hal_stamp,mirror,t = int(match[0]),match[1],int(match[2]),int(match[3])
      hent,hmean,hRMS = int(match[4]),float(match[5]),float(match[6])
      const,mean,sigma = float(match[7]),float(match[8]),float(match[9])
      if not mirror in stat_tuples:
        stat_tuples[mirror] = []
      stat_tuples[mirror].append((t,hent,hmean,hRMS,const,mean,sigma))
    return stat_tuples

  def get_flash_stats_tube(self):
    stat_tuples = {}
    ext_hs_regex = re.compile(r'(y\d{4}m\d{2}d\d{2}p\d{2})m(\d{2})t(\d{3})')
    for line in csv.reader(self.flash_stats_tube_file,delimiter=' '):
      i = int(line[0])
      groups = re.match(ext_hs_regex,line[1]).groups()
      hal_stamp,mirror,tube = groups[0],int(groups[1]),int(groups[2])
      t = int(line[2])
      hent,hmean,hRMS = int(line[3]),float(line[4]),float(line[5])
      const,mean,sigma = float(line[6]),float(line[7]),float(line[8])
      if not mirror in stat_tuples:
        stat_tuples[mirror] = {}
      if not tube in stat_tuples[mirror]:
        stat_tuples[mirror][tube] = {}
      stat_tuples[mirror][tube][hal_stamp] = (t,hent,hmean,hRMS,const,mean,sigma)
    return stat_tuples

  def plot_QDCB_tube_vs_time(self):
    flash_stats = self.get_flash_stats_tube()
    flash_hists = {}
    for mirror in flash_stats:
      if not mirror in flash_hists:
        flash_hists[mirror] = {}
      for tube in flash_stats[mirror]:
        if not tube in flash_hists[mirror]:
          flash_hists[mirror][tube] = TH2F('m%02dt%03d' % (mirror,tube),';time;m%02dt%03d QDCB LED mean' % (mirror,tube),*(self.time_bins + self.mean_bins))
        for entry in flash_stats[mirror][tube]:
          flash_hists[mirror][tube].Fill(entry[0]/1000. - self.t0,entry[2]) # hmean vs time
    self.write_plots('/home/findlay/data/plots/LED/QDCB_tube_vs_time_m06_t128-255.ps',flash_hists[6].values(),xaxis_time=True)

  def plot_LED_flash_stats(self):
    stats_hists = {}
    for line in self.flash_stats_cluster_file:
      match = re.match(self.flash_stats_cluster_regex,line).groups()
      i,hal_stamp,mirror,t = int(match[0]),match[1],int(match[2]),int(match[3])
      hent,hmean,hRMS = int(match[4]),float(match[5]),float(match[6])
      const,mean,sigma = float(match[7]),float(match[8]),float(match[9])
      if hent < 241 : continue
      time_diff = t/1000. - self.t0
      if not mirror in stats_hists:
        hent_hist = TH2F('mirror %d LED flash tubes' % mirror,';time [200[8|9]-mm-dd];QDCB tubes',*(self.time_bins + self.hent_bins))
        hmean_hist = TH2F('mirror %d LED flash hist mean' % mirror,';time [200[8|9]-mm-dd];QDCB hist mean',*(self.time_bins + self.mean_bins))
        hRMS_hist = TH2F('mirror %d LED flash hist RMS' % mirror,';time [200[8|9]-mm-dd];QDCB hist RMS',*(self.time_bins + self.RMS_bins))
        const_hist = TH2F('mirror %d LED flash constants' % mirror,';time [200[8|9]-mm-dd];QDCB normal constant',*(self.time_bins + self.const_bins))
        mean_hist = TH2F('mirror %d LED flash means' % mirror,';time [200[8|9]-mm-dd];QDCB normal mean',*(self.time_bins + self.mean_bins))
        sigma_hist = TH2F('mirror %d LED flash sigmas' % mirror,';time [200[8|9]-mm-dd];QDCB normal sigma',*(self.time_bins + self.RMS_bins))
        stats_hists[mirror] = {'hent':hent_hist,'hmean':hmean_hist,'hRMS':hRMS_hist,'const':const_hist,'mean':mean_hist,'sigma':sigma_hist}
      stats_hists[mirror]['hent'].Fill(time_diff,hent)
      stats_hists[mirror]['hmean'].Fill(time_diff,hmean)
      stats_hists[mirror]['hRMS'].Fill(time_diff,hRMS)
      stats_hists[mirror]['const'].Fill(time_diff,const)
      stats_hists[mirror]['mean'].Fill(time_diff,mean)
      stats_hists[mirror]['sigma'].Fill(time_diff,sigma)
    self.write_kind_plots(stats_hists,'/home/findlay/data/plots/LED/LED_flash_%s.ps',xaxis_time=True)

  def plot_LED_temps(self):
    LED_temp_hists = {}
    for hal_stamp in sorted(self.data['led355'].keys()):
      notice_branch = self.data['led355'][hal_stamp]['notice']
      notice = THPKT1_DST_NOTICE()
      notice_branch.SetAddress(AddressOf(notice))
      for entry in xrange(notice_branch.GetEntries()):
        notice_branch.GetEntry(entry)
        if notice.type == 10:
          try:
            tA,tB,tC,tD = [float(temp) for temp in re.match(self.LED_temp_regex,notice.text).groups()]
          except:
            continue
          mirror = notice.pktHdr_crate
          time_diff = convert_time(notice)/1000. - self.t0
          if not mirror in LED_temp_hists:
            hist_A = TH2F('mirror %d LED temp A' % mirror,';time [200[8|9]-mm-dd];temperature [K]',*(self.time_bins + self.LED_T_bins))
            hist_B = TH2F('mirror %d LED temp B' % mirror,';time [200[8|9]-mm-dd];temperature [K]',*(self.time_bins + self.LED_T_bins))
            hist_C = TH2F('mirror %d LED temp C' % mirror,';time [200[8|9]-mm-dd];temperature [K]',*(self.time_bins + self.LED_T_bins))
            hist_D = TH2F('mirror %d LED temp D' % mirror,';time [200[8|9]-mm-dd];temperature [K]',*(self.time_bins + self.T_bins))
            LED_temp_hists[mirror] = {'A':hist_A,'B':hist_B,'C':hist_C,'D':hist_D}
          LED_temp_hists[mirror]['A'].Fill(time_diff,tA)
          LED_temp_hists[mirror]['B'].Fill(time_diff,tB)
          LED_temp_hists[mirror]['C'].Fill(time_diff,tC)
          LED_temp_hists[mirror]['D'].Fill(time_diff,tD)
    self.write_kind_plots(LED_temp_hists,'/home/findlay/data/plots/LED/LED_temp_%s.ps',xaxis_time=True)

  def plot_temp_AB_average(self):
    average_hists = {}
    for hal_stamp in sorted(self.data['led355'].keys()):
      notice_branch = self.data['led355'][hal_stamp]['notice']
      notice = THPKT1_DST_NOTICE()
      notice_branch.SetAddress(AddressOf(notice))
      for entry in xrange(notice_branch.GetEntries()):
        notice_branch.GetEntry(entry)
        if notice.type == 10:
          try:
            tA,tB,tC,tD = [float(temp) for temp in re.match(self.LED_temp_regex,notice.text).groups()]
          except:
            continue
          mirror = notice.pktHdr_crate
          time_diff = convert_time(notice)/1000. - self.t0
          if not mirror in average_hists:
            average_hists[mirror] = TH2F('mirror %d temp A,B average' % mirror,';time [200[8|9]-mm-dd];temp (A + B)/2 [K]',*(self.time_bins + self.LED_T_bins))
          average_hists[mirror].Fill(time_diff,(tA + tB)/2.)
    self.write_plots('/home/findlay/data/plots/LED/LED_temp_AB_average.ps',average_hists.values(),xaxis_time=True)

  def plot_cluster_PTH(self):
    cluster_hists = {}
    for hal_stamp in sorted(self.data['led355'].keys()):
      notice_branch = self.data['led355'][hal_stamp]['notice']
      notice = THPKT1_DST_NOTICE()
      notice_branch.SetAddress(AddressOf(notice))
      for entry in xrange(notice_branch.GetEntries()):
        notice_branch.GetEntry(entry)
        if notice.type == 17:
          try:
            match = re.match(self.PTH_regex,notice.text).groups()
            mirror,press,temp,hum = int(match[0]),float(match[1]),float(match[2]),float(match[3])
          except:
            continue
          time_diff = convert_time(notice)/1000. - self.t0
          if not mirror in cluster_hists:
            P_hist = TH2F('mirror %d pressure' % mirror,';time [200[8|9]-mm-dd];cluster pressure [Pa]',*(self.time_bins + self.P_bins))
            T_hist = TH2F('mirror %d temperature' % mirror,';time [200[8|9]-mm-dd];cluster temperature [K]',*(self.time_bins + self.T_bins))
            H_hist = TH2F('mirror %d humidity' % mirror,';time [200[8|9]-mm-dd];cluster humidity [\%]',*(self.time_bins + self.H_bins))
            cluster_hists[mirror] = {'press':P_hist,'temp':T_hist,'hum':H_hist}
          cluster_hists[mirror]['press'].Fill(time_diff,press)
          cluster_hists[mirror]['temp'].Fill(time_diff,temp)
          cluster_hists[mirror]['hum'].Fill(time_diff,hum)
    self.write_kind_plots(cluster_hists,'/home/findlay/data/plots/LED/cluster_%s.ps',xaxis_time=True)

  def plot_AB_average_QDCB(self):
    temp_tuples = {} # LED temperatures by mirror
    temp_QDCB_hists = {}
    for hal_stamp in sorted(self.data['led355'].keys()):
      notice_branch = self.data['led355'][hal_stamp]['notice']
      notice = THPKT1_DST_NOTICE()
      notice_branch.SetAddress(AddressOf(notice))
      for entry in xrange(notice_branch.GetEntries()):
        notice_branch.GetEntry(entry)
        if notice.type == 10:
          try:
            tA,tB,tC,tD = [float(temp) for temp in re.match(self.LED_temp_regex,notice.text).groups()]
          except:
            continue
          if not notice.pktHdr_crate in temp_QDCB_hists:
            temp_QDCB_hists[notice.pktHdr_crate] = {}
            temp_tuples[notice.pktHdr_crate] = []
          temp_tuples[notice.pktHdr_crate].append((convert_time(notice),tA,tB,tC,tD))
    for mirror in temp_QDCB_hists:
      AB_mean_hist = TH2F('mirror %s temp A,B average hmean' % mirror,';temp (A + B)/2 [K];LED flash QDCB means',*(self.LED_AB_bins + self.mean_bins))
      AB_RMS_hist = TH2F('mirror %s temp A,B average hRMS' % mirror,';temp (A + B)/2 [K];LED flash QDCB RMSs',*(self.LED_AB_bins + self.RMS_bins))
      temp_QDCB_hists[mirror] = {'AB_average_QDCB_hmean':AB_mean_hist,'AB_average_QDCB_hRMS':AB_RMS_hist}
    stat_tuples = self.get_flash_stats_cluster()
    for mirror in stat_tuples:
      correlated_tuples = find_nearest_tuple(stat_tuples[mirror],temp_tuples[mirror])
      for tuple in correlated_tuples:
        # tuple - 0:t 1:hent 2:hmean 3:hRMS 4:const 5:mean 6:sigma 7:t 8:tA 9:tB 10:tC 11:tD
        temp_QDCB_hists[mirror]['AB_average_QDCB_hmean'].Fill((tuple[8] + tuple[9])/2.,tuple[2]) # (tA + tB)/2,hmean
        temp_QDCB_hists[mirror]['AB_average_QDCB_hRMS'].Fill((tuple[8] + tuple[9])/2.,tuple[3]) # (tA + tB)/2,hRMS
    self.write_kind_plots(temp_QDCB_hists,'/home/findlay/data/plots/LED/%s.ps')

  def plot_temp_QDCB(self):
    temp_tuples = {} # LED temperatures by mirror
    temp_QDCB_hists = {}
    for hal_stamp in sorted(self.data['led355'].keys()):
      notice_branch = self.data['led355'][hal_stamp]['notice']
      notice = THPKT1_DST_NOTICE()
      notice_branch.SetAddress(AddressOf(notice))
      for entry in xrange(notice_branch.GetEntries()):
        notice_branch.GetEntry(entry)
        if notice.type == 10:
          try:
            tA,tB,tC,tD = [float(temp) for temp in re.match(self.LED_temp_regex,notice.text).groups()]
          except:
            continue
          if not notice.pktHdr_crate in temp_QDCB_hists:
            temp_QDCB_hists[notice.pktHdr_crate] = {}
            temp_tuples[notice.pktHdr_crate] = []
          temp_tuples[notice.pktHdr_crate].append((convert_time(notice),tA,tB,tC,tD))
    for mirror in temp_QDCB_hists:
      A_mean_hist = TH2F('mirror %s LED temp A hmean' % mirror,';LED temperature A [K];LED flash QDCB means',*(self.LED_T_bins + self.mean_bins))
      B_mean_hist = TH2F('mirror %s LED temp B hmean' % mirror,';LED temperature B [K];LED flash QDCB means',*(self.LED_T_bins + self.mean_bins))
      C_mean_hist = TH2F('mirror %s LED temp C hmean' % mirror,';LED temperature C [K];LED flash QDCB means',*(self.LED_T_bins + self.mean_bins))
      D_mean_hist = TH2F('mirror %s LED temp D hmean' % mirror,';LED temperature D [K];LED flash QDCB means',*(self.T_bins + self.mean_bins))
      A_RMS_hist = TH2F('mirror %s LED temp A hRMS' % mirror,';LED temperature A [K];LED flash QDCB RMSs',*(self.LED_T_bins + self.RMS_bins))
      B_RMS_hist = TH2F('mirror %s LED temp B hRMS' % mirror,';LED temperature B [K];LED flash QDCB RMSs',*(self.LED_T_bins + self.RMS_bins))
      C_RMS_hist = TH2F('mirror %s LED temp C hRMS' % mirror,';LED temperature C [K];LED flash QDCB RMSs',*(self.LED_T_bins + self.RMS_bins))
      D_RMS_hist = TH2F('mirror %s LED temp D hRMS' % mirror,';LED temperature D [K];LED flash QDCB RMSs',*(self.T_bins + self.RMS_bins))
      temp_QDCB_hists[mirror] = {'temp_A_QDCB_hmean':A_mean_hist,'temp_B_QDCB_hmean':B_mean_hist,
          'temp_C_QDCB_hmean':C_mean_hist,'temp_D_QDCB_hmean':D_mean_hist,
          'temp_A_QDCB_hRMS':A_RMS_hist,'temp_B_QDCB_hRMS':B_RMS_hist,
          'temp_C_QDCB_hRMS':C_RMS_hist,'temp_D_QDCB_hRMS':D_RMS_hist}
    stat_tuples = self.get_flash_stats_cluster()
    for mirror in stat_tuples:
      correlated_tuples = find_nearest_tuple(stat_tuples[mirror],temp_tuples[mirror])
      for tuple in correlated_tuples:
        # tuple - 0:t 1:hent 2:hmean 3:hRMS 4:const 5:mean 6:sigma 7:t 8:tA 9:tB 10:tC 11:tD
        temp_QDCB_hists[mirror]['temp_A_QDCB_hmean'].Fill(tuple[8],tuple[2]) # tA,hmean
        temp_QDCB_hists[mirror]['temp_B_QDCB_hmean'].Fill(tuple[9],tuple[2]) # tB,hmean
        temp_QDCB_hists[mirror]['temp_C_QDCB_hmean'].Fill(tuple[10],tuple[2]) # tC,hmean
        temp_QDCB_hists[mirror]['temp_D_QDCB_hmean'].Fill(tuple[11],tuple[2]) # tD,hmean
        temp_QDCB_hists[mirror]['temp_A_QDCB_hRMS'].Fill(tuple[8],tuple[3]) # tA,hRMS
        temp_QDCB_hists[mirror]['temp_B_QDCB_hRMS'].Fill(tuple[9],tuple[3]) # tB,hRMS
        temp_QDCB_hists[mirror]['temp_C_QDCB_hRMS'].Fill(tuple[10],tuple[3]) # tC,hRMS
        temp_QDCB_hists[mirror]['temp_D_QDCB_hRMS'].Fill(tuple[11],tuple[3]) # tD,hRMS
    self.write_kind_plots(temp_QDCB_hists,'/home/findlay/data/plots/LED/%s.ps')

  def plot_PTH_QDCB(self):
    PTH_tuples = {} # cluster PTH by mirror
    PTH_QDCB_hists = {}
    for hal_stamp in sorted(self.data['led355'].keys()):
      notice_branch = self.data['led355'][hal_stamp]['notice']
      notice = THPKT1_DST_NOTICE()
      notice_branch.SetAddress(AddressOf(notice))
      for entry in xrange(notice_branch.GetEntries()):
        notice_branch.GetEntry(entry)
        if notice.type == 17:
          match = re.match(self.PTH_regex,notice.text)
          if match == None : continue
          groups = match.groups()
          mirror,press,temp,hum = int(groups[0]),float(groups[1]),float(groups[2]),float(groups[3])
          if not mirror in PTH_QDCB_hists:
            PTH_QDCB_hists[mirror] = {}
            PTH_tuples[mirror] = []
          PTH_tuples[mirror].append((convert_time(notice),press,temp,hum))
    for mirror in PTH_QDCB_hists:
      P_mean_hist = TH2F('mirror %s cluster press hmean' % mirror,';cluster pressure [Pa];LED flash QDCB means',*(self.P_bins + self.mean_bins))
      T_mean_hist = TH2F('mirror %s cluster temp hmean' % mirror,';cluster temperature [K];LED flash QDCB means',*(self.T_bins + self.mean_bins))
      H_mean_hist = TH2F('mirror %s cluster hum hmean' % mirror,';cluster humidity [\%];LED flash QDCB means',*(self.H_bins + self.mean_bins))
      P_RMS_hist = TH2F('mirror %s cluster press hRMS' % mirror,';cluster pressure [Pa];LED flash QDCB RMSs',*(self.P_bins + self.RMS_bins))
      T_RMS_hist = TH2F('mirror %s cluster temp hRMS' % mirror,';cluster temperature [K];LED flash QDCB RMSs',*(self.T_bins + self.RMS_bins))
      H_RMS_hist = TH2F('mirror %s cluster hum hRMS' % mirror,';cluster humidity [\%];LED flash QDCB RMSs',*(self.H_bins + self.RMS_bins))
      PTH_QDCB_hists[mirror] = {'cluster_press_QDCB_hmean':P_mean_hist,'cluster_temp_QDCB_hmean':T_mean_hist,
          'cluster_hum_QDCB_hmean':H_mean_hist,'cluster_press_QDCB_hRMS':P_RMS_hist,
          'cluster_temp_QDCB_hRMS':T_RMS_hist,'cluster_hum_QDCB_hRMS':H_RMS_hist}
    stat_tuples = self.get_flash_stats_cluster()
    for mirror in stat_tuples:
      correlated_tuples = find_nearest_tuple(stat_tuples[mirror],PTH_tuples[mirror])
      for tuple in correlated_tuples:
        # tuple 0:t 1:hent 2:hmean 3:hRMS 4:const 5:mean 6:sigma 7:t 8:press 9:temp 10:hum
        PTH_QDCB_hists[mirror]['cluster_press_QDCB_hmean'].Fill(tuple[8],tuple[2]) # press,hmean
        PTH_QDCB_hists[mirror]['cluster_temp_QDCB_hmean'].Fill(tuple[9],tuple[2]) # temp,hmean
        PTH_QDCB_hists[mirror]['cluster_hum_QDCB_hmean'].Fill(tuple[10],tuple[2]) # hum,hmean
        PTH_QDCB_hists[mirror]['cluster_press_QDCB_hRMS'].Fill(tuple[8],tuple[3]) # press,hRMS
        PTH_QDCB_hists[mirror]['cluster_temp_QDCB_hRMS'].Fill(tuple[9],tuple[3]) # temp,hRMS
        PTH_QDCB_hists[mirror]['cluster_hum_QDCB_hRMS'].Fill(tuple[10],tuple[3]) # hum,hRMS
    self.write_kind_plots(PTH_QDCB_hists,'/home/findlay/data/plots/LED/%s.ps')

class CalibHV_QDCB(LEDEnv):
  '''calibrate MDFD'''

  def plot_m06_tube_QDCB_means_subcl(self):
    stat_tuples = self.get_flash_stats_tube()
    mirror = 6
    mean_QDCB_hists = {mirror:{}}
    for tube in stat_tuples[mirror]:
      subcl = tube/16 + 1
      if not subcl in mean_QDCB_hists[mirror]:
        subcl_hist = TProfile('m%02dsubcl%02d' % (mirror,subcl),'m06 mean QDCBs;tube number;tube QDCB mean [counts]',16,16*(subcl - 1),16*(subcl - 1) + 16)
        subcl_hist.SetFillColor(subcl)
        subcl_hist.SetLineColor(subcl)
        mean_QDCB_hists[mirror][subcl] = subcl_hist
      for hal_stamp in stat_tuples[mirror][tube]:
        mean_QDCB_hists[mirror][subcl].Fill(tube,stat_tuples[mirror][tube][hal_stamp][2])
    self.write_plots('/home/findlay/data/plots/QDCB/m06_tube_QDCB_means_subcl.ps',mean_QDCB_hists[mirror].values())

  def plot_m06_tube_HV_means_subcl(self):
    tube_HV_mean_hists = {}
    for hal_stamp in sorted(self.data['led355'].keys()):
      volts_branch = self.data['led355'][hal_stamp]['volts']
      volts = THPKT1_DST_VOLTS()
      volts_branch.SetAddress(AddressOf(volts))
      for entry in xrange(volts_branch.GetEntries()):
        volts_branch.GetEntry(entry)
        mirror = volts.pktHdr_crate
        if not mirror == 6 : continue
        if not mirror in tube_HV_mean_hists:
          tube_HV_mean_hists[mirror] = {}
        for tube in xrange(len(volts.hv)):
          subcl = tube/16 + 1
          if not subcl in tube_HV_mean_hists[mirror]:
            subcl_hist = TProfile('m%02dsubcl%02d' % (mirror,subcl),'m06 mean QDCBs;tube number;tube QDCB mean [counts]',16,16*(subcl - 1),16*(subcl - 1) + 16)
            subcl_hist.SetFillColor(subcl)
            subcl_hist.SetLineColor(subcl)
            tube_HV_mean_hists[mirror][subcl] = subcl_hist
          tube_HV_mean_hists[mirror][subcl].Fill(tube,volts.hv[tube])
    self.write_plots('/home/findlay/data/plots/HV/m06_tube_HV_means_subcl.ps',tube_HV_mean_hists[6].values())

  def plot_m06_tube_QDCB_means_mirror(self):
    stat_tuples = self.get_flash_stats_tube()
    mean_QDCB_stack = THStack('m06 mean QDCBs','m06 mean QDCBs;tube number;tube QDCB mean [V]')
    hist_legend = TLegend(0.89,0.89,0.90,0.90)
    for mirror in stat_tuples:
      if not mirror == 6 : continue
      for tube in stat_tuples[mirror]:
        subcl = tube/16 + 1
        subcl_hist = TProfile('m%02dt%03d' % (mirror,tube),'m06 mean QDCBs;tube number;tube QDCB mean [V]',self.tubes + 1,0,self.tubes)
        subcl_hist.SetFillColor(subcl)
        subcl_hist.SetLineColor(subcl)
        for hal_stamp in stat_tuples[mirror][tube]:
          subcl_hist.Fill(tube,stat_tuples[mirror][tube][hal_stamp][2])
#        hist_legend.AddEntry(subcl_hist,'subcl %d' % subcl,'f')
        mean_QDCB_stack.Add(subcl_hist)
    self.write_plot('/home/findlay/data/plots/QDCB/m06_tube_QDCB_means.ps',mean_QDCB_stack,legend=hist_legend)
    self.convert_plot_file('/home/findlay/data/plots/QDCB/m06_tube_QDCB_means.ps')

  def plot_m06_tube_HV_means_mirror(self):
    tube_HV_mean_hists = {}
    mean_HV_stack = THStack('m06 mean HVs','m06 mean HVs; ; ')
    hist_legend = TLegend(0.89,0.89,0.90,0.90)
    for hal_stamp in sorted(self.data['led355'].keys()):
      volts_branch = self.data['led355'][hal_stamp]['volts']
      volts = THPKT1_DST_VOLTS()
      volts_branch.SetAddress(AddressOf(volts))
      for entry in xrange(volts_branch.GetEntries()):
        volts_branch.GetEntry(entry)
        mirror = volts.pktHdr_crate
        if not mirror == 6 : continue
        if not mirror in tube_HV_mean_hists:
          tube_HV_mean_hists[mirror] = {}
        for tube in xrange(len(volts.hv)):
          subcl = tube/16 + 1
          if not subcl in tube_HV_mean_hists[mirror]:
            subcl_hist = TProfile('m%02dt%03d' % (mirror,tube),';tube number;tube HV mean [V]',self.tubes + 1,0,self.tubes)
            subcl_hist.SetFillColor(subcl)
            subcl_hist.SetLineColor(subcl)
            mean_HV_stack.Add(subcl_hist)
#            hist_legend.AddEntry(subcl_hist,'subcl %d' % subcl,'f')
            tube_HV_mean_hists[mirror][subcl] = subcl_hist
          tube_HV_mean_hists[mirror][subcl].Fill(tube,volts.hv[tube])
    self.write_plot('/home/findlay/data/plots/HV/m06_tube_HV_means.ps',mean_HV_stack,legend=hist_legend)
    self.convert_plot_file('/home/findlay/data/plots/HV/m06_tube_HV_means.ps')

  def plot_HV_vs_temp(self):
    pth = {}
    i = 0
    for hal_stamp in sorted(self.data['led355'].keys()):
      pth[hal_stamp] = {}
      start = self.get_run_start(hal_stamp,'led355')
      notice_branch = self.data['led355'][hal_stamp]['notice']
      notice = THPKT1_DST_NOTICE()
      notice_branch.SetAddress(AddressOf(notice))
      for entry in xrange(notice_branch.GetEntries()):
        notice_branch.GetEntry(entry)
        if notice.type == 17:
          match = re.match(self.PTH_regex,notice.text)
          if match == None : continue
          groups = match.groups()
          mirror,press,temp,hum = int(groups[0]),float(groups[1]),float(groups[2]),float(groups[3])
          if not mirror in pth[hal_stamp] :
            pth[hal_stamp][mirror] = {'press':[],'temp':[],'hum':[]}
          pth[hal_stamp][mirror]['press'].append(press)
          pth[hal_stamp][mirror]['temp'].append(temp)
          pth[hal_stamp][mirror]['hum'].append(hum)
    HV_vs_temp_hists = {}
    for hal_stamp in sorted(self.data['led355'].keys()):
      start = self.get_run_start(hal_stamp,'led355')
      volts_branch = self.data['led355'][hal_stamp]['notice']
      volts = THPKT1_DST_VOLTS()
      volts_branch.SetAddress(AddressOf(volts))
      for entry in xrange(volts_branch.GetEntries()):
        volts_branch.GetEntry(entry)
        mirror = volts.pktHdr_crate
        if not mirror == 6 : continue
        if not mirror in HV_vs_temp_hists:
          HV_vs_temp_hists[mirror] = {}
        for tube in xrange(len(volts.hv)):
          if tube > 128 : continue
          #print hal_stamp,mirror,tube,volts.hv[tube] ; raw_input()
          if not tube in HV_vs_temp_hists[mirror]:
            HV_vs_temp_hists[mirror][tube] = TH2I('m%02dt%03d' % (mirror,tube),';temperature [K];HV [V]',*(self.T_bins + (512,0,20)))
          HV_vs_temp_hists[mirror][tube].Fill(mean(pth[hal_stamp][mirror]['temp']),volts.hv[tube])
    self.write_plots('/home/findlay/data/plots/LED/HV_vs_temp_m06_t000-127.ps',HV_vs_temp_hists[6].values())

  def plot_QDCB_vs_HV(self):
    stat_tuples = self.get_flash_stats_tube()
    QDCB_vs_HV_hists = {}
    for hal_stamp in sorted(self.data['led355'].keys()):
      start = self.get_run_start(hal_stamp,'led355')
      volts_branch = self.data['led355'][hal_stamp]['volts']
      volts = THPKT1_DST_VOLTS()
      volts_branch.SetAddress(AddressOf(volts))
      for entry in xrange(volts_branch.GetEntries()):
        volts_branch.GetEntry(entry)
        mirror = volts.pktHdr_crate
        if not mirror == 6 : continue
        if not mirror in QDCB_vs_HV_hists:
          QDCB_vs_HV_hists[mirror] = {}
        for tube in xrange(len(volts.hv)):
          if tube > 128 : continue
          if not tube in QDCB_vs_HV_hists[mirror]:
            QDCB_vs_HV_hists[mirror][tube] = TH2I('m%02dt%03d' % (mirror,tube),';HV [V];QDCB LED means',*(self.volts_bins + self.mean_bins))
          if hal_stamp in stat_tuples[mirror][tube]:
            QDCB_vs_HV_hists[mirror][tube].Fill(volts.hv[tube],stat_tuples[mirror][tube][hal_stamp][2]) # HV,QDCB mean
    self.write_plots('/home/findlay/data/plots/LED/QDCB_tube_vs_HV_m06_t000-127.ps',QDCB_vs_HV_hists[6].values())

class CalibTDC_QDCB():
  '''try to model all of the different noises seen in LED flashes so that we can reliably subtract them'''

  class ClusterModel(Plot):
    '''
    Assume that low brightness flashes trigger the TDC so that it either:
     - includes some extra noise (pedestal noise? junk LED flash?),
     - triggers early and cuts off part of the LED light,
     - triggers early and cuts off all of the LED light.
    Generate a pedestal distribution for each tube from the entire noise-closed set and use these distributions in the convolution of the cumulative distribution of QDC data when the data is asymmetrically longer below the mean.

    When HAL collects data in the course of an 'event' it reduces the three dimensions of analog information (tube,time,signal) into a packet of snapshot information:
     - a 2D projection of the QDC distribution across all tubes participating in the event in two scales (QDCA,QDCB),
     - a scalar indicating the duration of the event,
     - an accurate timepoint indicating the beginning of the event.
    This is the information we have to work with in order to verify our model.

    Need to estimate N_pe = k(mu/sigma)^2, not mu and sigma individually, but mu/sigma, (mu/sigma)^2 if possible from statistical data.  This provides one more degree of freedom
    (mu/sigma)^2 = 1/[M(sum(x*x)/sum^2(x) - 2) + 1]
    '''
    '''First, lets analyze the general trends in QDC and TDC dists'''
    '''
    provisional calibration sketch:  For each led355 part for each cluster
    - fit cluster TDC: cut on 2 sigma
    - fit cluster QDCB
    - remove tube data whose means are without 2 sigma of the cluster TDC,QDCB dists
    - refit cluster TDC,QDCB dists
    - extrapolate N_pe from mu,sigma
    '''
    def observe_flashes(self):
      self.canvas.cd()
      self.canvas.Divide(2,2)
      batch = True
      if '-d' in sys.argv : batch = False
      for hal_stamp in sorted(self.data['led355'].keys()):
        print
        event_branch = self.data['led355'][hal_stamp]['event']
        event = THPKT1_DST_EVENT()
        event_branch.SetAddress(AddressOf(event))
        for entry in xrange(event_branch.GetEntries()):
          event_branch.GetEntry(entry)
          if mean(event.qdcB) < 800:
            if batch:
              mirror = event.pktHdr_crate
              info_tuple = (hal_stamp,mirror,entry,
                  max(event.qdcB),mean(event.qdcB),std(event.qdcB),min(event.qdcB),
                  max(event.tdc),mean(event.tdc),std(event.tdc),min(event.tdc))
              print '%s mirror %02d event %04d  QDCB(max %4d mean %4d var %4d min %4d)  TDC(max %4d mean %4d var %4d min %4d)' % info_tuple
              last_hal_stamp,last_mirror,last_entry = hal_stamp,mirror,entry
              if hal_stamp == last_hal_stamp and mirror < last_mirror and entry > last_entry:
                print '\n  mirror %02d event %04d out of order\n'
            else:
              h1 = TH2F('%sm%02dentry%da' % (hal_stamp,mirror,entry),';m%02d tube;QDCB' % mirror,256,0,255,512,0,2048)
              h2 = TH2F('%sm%02dentry%db' % (hal_stamp,mirror,entry),';TDC;QDCB',1024,0,4352,512,0,2048)
              h3 = TH2F('%sm%02dentry%dc' % (hal_stamp,mirror,entry),';m%02d tube;TDC' % mirror,256,0,256,1024,0,4352)
              for i in xrange(event.ntubes):
                h1.Fill(event.tube_num[i],event.qdcB[i])
                h2.Fill(event.tdc[i],event.qdcB[i])
                h3.Fill(event.tube_num[i],event.tdc[i])
              print hal_stamp,mirror,entry
              self.display_histogram(h1,self.canvas.cd(1),pause=False)
              self.display_histogram(h2,self.canvas.cd(2),pause=False)
              self.display_histogram(h3,self.canvas.cd(3),pause=False)
              raw_input()

    def observe_parts(self):
      self.canvas.cd()
      self.canvas.Divide(1,2)
      batch = True
      if '-d' in sys.argv : batch = False
      part_data = []
      part_hists = []
      for hal_stamp in sorted(self.data['led355'].keys()):
        if hal_stamp > HalStamp('y2009m07d01p01'):
          event_branch = self.data['led355'][hal_stamp]['event']
          event = THPKT1_DST_EVENT()
          event_branch.SetAddress(AddressOf(event))
          part_data.append({'q':{},'t':{}})
          part_hists.append({'hq':{},'ht':{}})
          if batch:
            q = part_data[-1]['q']
            t = part_data[-1]['t']
            for i in xrange(1,self.mirrors + 1):
              q[i] = []
              t[i] = []
            for entry in xrange(event_branch.GetEntries()):
              event_branch.GetEntry(entry)
              mirror = event.pktHdr_crate
              for i in xrange(event.ntubes):
                q[mirror].append(event.qdcB[i])
                t[mirror].append(event.tdc[i])
            for mirror in xrange(1,self.mirrors + 1):
              if len(q[mirror]) != 0:
                print self.find_peak(q[mirror])
                raw_input()
          else:
            hq = part_hists[-1]['hq'] ; cq = self.canvas.cd(1)
            ht = part_hists[-1]['ht'] ; ct = self.canvas.cd(2) ; ct.SetLogy()
            for m in xrange(1,self.mirrors + 1):
              hq[m] = TH1F('%sm%02dq' % (hal_stamp,m),';QDCB;counts',128,0,4096)
              ht[m] = TH1F('%sm%02dt' % (hal_stamp,m),';TDC;counts',128,0,4096)
            for entry in xrange(event_branch.GetEntries()):
              event_branch.GetEntry(entry)
              mirror = event.pktHdr_crate
              for i in xrange(event.ntubes):
                hq[mirror].Fill(event.qdcB[i] + random())
                ht[mirror].Fill(event.tdc[i] + random())
            print hal_stamp
            for mirror in xrange(1,self.mirrors + 1):
              if hq[mirror].GetEntries() != 0:
              	cq.cd() ; hq[mirror].Fit('gaus','M Q') ; self.display_histogram(hq[mirror],cq,pause=False)
                ct.cd() ; ht[mirror].Fit('gaus','M Q') ; self.display_histogram(ht[mirror],ct,pause=False)
                raw_input()

    def provisional_calibrate(self):
      for hal_stamp in sorted(self.data['led355'].keys()):
        if hal_stamp > HalStamp('y2009m07d01p01'):
          event_branch = self.data['led355'][hal_stamp]['event']
          event = THPKT1_DST_EVENT()
          event_branch.SetAddress(AddressOf(event))
          part_data = {}
          for m in xrange(1,self.mirrors + 1):
            part_data[m] = {}
            for t in xrange(self.tubes):
              part_data[m][t] = []
          for entry in xrange(event_branch.GetEntries()):
            event_branch.GetEntry(entry)
            mirror = event.pktHdr_crate
            for t in xrange(event.ntubes):
              part_data[mirror][event.tube_num[t]].append([True,event.tdc[t],event.qdcB[t]]) # keep,tdc,qdcb

          '''fit cluster TDCs, QDCBs, and tube QDCBs, cutting events as appropriate'''
          self.fit_cluster_TDCs(hal_stamp,part_data)
          QDCB_cluster_fit_data = self.fit_cluster_QDCBs(hal_stamp,part_data)
          QDCB_tube_fit_data = self.fit_tube_QDCBs(hal_stamp,part_data,QDCB_cluster_fit_data)
          self.tally(part_data)

          '''refit cluster QDCB and tube QDCB'''
          QDCB_cluster_fit_data = self.fit_cluster_QDCBs(hal_stamp,part_data)
          QDCB_tube_fit_data = self.fit_tube_QDCBs(hal_stamp,part_data,QDCB_cluster_fit_data)
          self.tally(part_data)

    def fit_cluster_TDCs(self,hal_stamp,part_data):
      for m in part_data:
        h = TH1F('%sm%02dt' % (hal_stamp,m),';TDC;counts',128,0,4096)
        for t in part_data[m]:
          for i in part_data[m][t]:
            if i[0] == True:
              h.Fill(i[1])
        h.Fit('gaus','L M Q')
        gaus = h.GetListOfFunctions().FindObject('gaus')
        const,mean,var = gaus.GetParameter(0),gaus.GetParameter(1),gaus.GetParameter(2)
        for t in part_data[m]:
          for i in part_data[m][t]:
            if fabs(i[1] - mean) > 4*var:
              i[0] = False
        #self.display_histogram(h)

    def fit_cluster_QDCBs(self,hal_stamp,part_data):
      fit_data = {}
      for m in part_data:
        h = TH1F('%sm%02dq' % (hal_stamp,m),';QDCB;counts',128,0,4096)
        for t in part_data[m]:
          for i in part_data[m][t]:
            if i[0] == True:
              h.Fill(i[2])
        h.Fit('gaus','L M Q')
        gaus = h.GetListOfFunctions().FindObject('gaus')
        fit_data[m] = (gaus.GetParameter(0),gaus.GetParameter(1),gaus.GetParameter(2))
        #self.display_histogram(h)
      return fit_data

    def fit_tube_QDCBs(self,hal_stamp,part_data,cluster_data):
      fit_data = {}
      for m in part_data:
        fit_data[m] = {}
        for t in part_data[m]:
          h = TH1F('%sm%02dt%03dq' % (hal_stamp,m,t),';QDCB;counts',128,0,4096)
          # fill hist
          for i in part_data[m][t]:
            if i[0] == True:
              h.Fill(i[2])
          # make cuts
          if h.GetEntries() != 0:
            status = h.Fit('gaus','L M Q')
            gaus = h.GetListOfFunctions().FindObject('gaus')
            fit_data[m][t] = (gaus.GetParameter(0),gaus.GetParameter(1),gaus.GetParameter(2))
            if fabs(fit_data[m][t][1] - cluster_data[m][1]) > 3*cluster_data[m][2]: # |tube_mean - cluster_mean| > 3*cluster_var
              for i in part_data[m][t]:
                i[0] = False
      return fit_data

    def tally(self,part_data):
      num = {}
      for m in xrange(1,self.mirrors + 1):
        num[m] = 0
        for t in part_data[m]:
          for i in part_data[m][t]:
            if i[0] == True:
              num[m] += 1
        print m,num[m]
      return num
      
    def find_peak(self,array):
      histogram = {}
      for bin in xrange(0,4096):
        histogram[bin] = 0
      for i in array:
        histogram[i] += 1
      print sorted(histogram.values()) ; raw_input()

    def generate_data(self):
      pass

def main():
  cm = CalibTDC_QDCB.ClusterModel()
  cm.provisional_calibrate()

if __name__ == '__main__' : main()
