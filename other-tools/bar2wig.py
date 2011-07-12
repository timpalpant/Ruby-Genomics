#!/usr/bin/python

import sys, os
from optparse import OptionParser
from MAT import AffyFileParser
def parse(fname, level, step):
    
    if not os.path.isfile(fname):
        print 'Error: could not open %s' % fname
        sys.exit('Stopping')

    oname = fname + '.wig'
    
    fp = open(oname, 'wt')
    fp.write( 'track type=wiggle_0\n' )

    # mat file parsing parsing
    from MAT import AffyFileParser
    bar  =  AffyFileParser.CBARFileWriter()
    bar.SetFileName(fname)
    if not bar.Read():
        raise Exception('Unable to properly read %s' % fname )

    seq  = AffyFileParser.CGDACSequenceResultItem()
    data = AffyFileParser.BarSequenceResultData()

    seq_range = range(bar.GetNumberSequences() - 1)
    
    for s in seq_range:
        bar.GetResults(s, seq)
        name  = seq.GetName()        
        for i in range( 0, seq.GetNumberDataPoints()-1, step ):
            seq.GetData(i, 0, data)
            x = data.iValue
            seq.GetData(i, 1, data)
            y = data.fValue
            if y >= level:                
                out ='%s\t%s\t%s\t%s\n' % ( name, x-1, x, y )
                fp.write( out )

    print '...saved bar data to: %s' % oname
                    

def usage():
    "Prints usage information"
    
    help = """
    
    USAGE: 
    
    bar2wig -t threshold -s step input_file

    PARAMETERS:
    
        - threshold is the minimal signal strength (default value 0.0)
        - step specifies the increment of the indices when sampling the  
          data and must be an integer >= 1 (default 1)
        
    Unspecified parameters will take the default values.

    EXAMPLE:
        
        bar2wig -t 1.5 -s 5 somedata.bar

    Will create the somedata.bar.wig file with threshold of 1.5
    while examining every 5th measurement. 

    SEE ALSO: wig2peak, wigalign
    """
   
    return help

if __name__ == '__main__':
    #sys.argv.extend( [ 'bar2text.py', '2' ] )
    try:
        import psyco
        psyco.full()
    except:
        pass

    parser = OptionParser()
    parser.set_usage( usage() )

    parser.add_option("-t", type="float", dest="level", default=0 )
    parser.add_option("-s", type="int", dest="step", default=1)
    (options, args) = parser.parse_args()

    if len(args)!=1:
        parser.print_usage()
        sys.exit(-1)
    else:
        fname = args[-1]

    if options.step < 1:
        parser.print_usage()
        print 'Incorrect step option %s. It must be greater or equal to 1' % options.step
        sys.exit()

    print '\n...executing bar2wig with threshold=%s, step=%d' % (options.level, options.step) 

    parse(fname=fname, level=options.level, step=options.step)

    print '...done\n'
