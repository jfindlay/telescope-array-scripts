#!/usr/bin/env /usr/bin/python
import os,sys,re,time
from datetime import date,datetime
from subprocess import call,Popen,PIPE

from calib import LEDFlashes

class MDData:
  def __init__(self):
    self.src_dir = '/tmp/middle_drum/'
    self.dest_dir = '/home/findlay/data/middle_drum/'
    self.hal_files = []
    self.led355_files = {}
    self.noise_closed_files = {}
    self.hvcalib_files = {}

  def mount_server(self):
    if not os.path.isdir(self.src_dir):
      os.makedirs(self.src_dir)
    if not os.path.ismount(self.src_dir):
      call(['sshfs','tamember@tadserv.physics.utah.edu:/tadserv2_4/tamd/',
          self.src_dir])

  def unmount_server(self):
    if os.path.ismount(self.src_dir):
      call(['fusermount','-u',self.src_dir])
    if os.path.isdir(self.src_dir):
      os.rmdir(self.src_dir)

  def collect_files(self):
    find_command = 'find %s -type f | grep -E "[0-9]{8}"' % self.src_dir
    find = Popen(find_command,stdout=PIPE,shell=True)
    self.hal_files = re.split('\n',find.stdout.read())[:-1]
    find.stdout.close()
    for name in self.hal_files:
      try:
        hal_stamp = re.search('(y\d{4}m\d{2}d\d{2}p\d{2})',name).group(0)
        led355,noise_closed,hvcalib = re.search('(led355.hal)|(noise-closed.hal)|(hvcalib.hal)$',name).groups()
      except:
        continue
      if led355:
        self.led355_files[hal_stamp] = name
      if noise_closed:
        self.noise_closed_files[hal_stamp] = name
      elif hvcalib:
        self.hvcalib_files[hal_stamp] = name

  # TODO: use method decorator
  def convert_files(self):
    for hal_stamp in sorted(self.led355_files.keys()):
      self.convert_file(self.led355_files[hal_stamp])
    for hal_stamp in sorted(self.noise_closed_files.keys()):
      self.convert_file(self.noise_closed_files[hal_stamp])
    for hal_stamp in sorted(self.hvcalib_files.keys()):
      self.convert_file(self.hvcalib_files[hal_stamp])

  def convert_file(self,hal_file):
    path,name = os.path.split(re.split(self.src_dir,hal_file)[1])
    base = os.path.splitext(name)[0]
    dest = os.path.join(self.dest_dir,path)
    if not os.path.isdir(dest):
      os.makedirs(dest)
    root_file = '%s.root' % os.path.join(dest,base)
    if not os.path.isfile(root_file):
      call(('/home/findlay/bin/hal2root',hal_file,root_file))
      print hal_file

def main():
  md = MDData()
  md.mount_server()
  md.collect_files()
  md.convert_files()
  md.unmount_server()

  lf = LEDFlashes()
  lf.make_flash_stats('/tmp/flash_stats_cluster.gz','/tmp/flash_stats_tube.gz')

if __name__ == '__main__' : main()
