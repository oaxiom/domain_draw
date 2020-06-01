
import matplotlib.cm as cm
import numpy

def collate_family_defining(filename):
    """
    scan a filename and list all 'FAMILY-DEFINING' features
    """
    
    oh = open(filename, "rU")
    
    doms = []
    for line in oh:
        if "FAMILY-DEFINING" in line:
            line = line.strip().split()
            doms.append(line[5])
    
    set_of_doms = set(doms)
    for d in set_of_doms:
        print("%s\t%s" % (doms.count(d), d))
    
    cols = cm.Paired(numpy.arange(len(set_of_doms))/ (len(set_of_doms)*1.0))[:,:-1]
    print(cols)
    
    print(cols)
    
    print()
    print("Suggested ColourMap:")
    print("col_map = {")
    for i, d in enumerate(set_of_doms):
        print("\t'%s': '%s'," % (d, '#%02x%02x%02x' % tuple(cols[i]*255)))
    print("\t}")
    
    oh.close()