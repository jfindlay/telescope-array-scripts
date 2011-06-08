import os,sys,re,csv

stan_reader = csv.reader(open('/home/findlay/data/pegs/newpegs.dat'),delimiter=' ')
tyler_reader = csv.reader(open('/home/findlay/data/pegs/newpegs.tyler.dat'),delimiter=' ')

def convert(reader=tyler_reader,chart_file='/home/findlay/data/pegs/tyler.pegs'):
  pegs = {}
  for row in reader:
    mirror,tube,peg = [int(i) for i in row]
    tube -= 1
    if not tube in pegs:
      pegs[tube] = {}
    pegs[tube][mirror] = peg
  
  f = file(chart_file,'w')
  f.write('Tube\tSub-c #\tTube #\tM1\tM2\tM3\tM4\tM5\tM6\tM7\tM8\tM9\tM10\tM11\tM12\tM13\tM14\n')
  for tube in pegs:
    line = '%d\t%d\t%d\t' % (tube,tube/16 + 1,tube%16 + 1)
    for mirror in pegs[tube]:
      if mirror != 14:
        line += '%d\t' % pegs[tube][mirror]
      else:
        line += '%d\n' % pegs[tube][mirror]
    f.write(line)

if __name__ == '__main__' : convert()
