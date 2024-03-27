import os
import sys
import pandas as pd

'''
HLA LOH interpretation criteria
1) A copy number < 0.5, is classified as subject to loss, and thereby indicative of LOH.
2) Allelic imbalance is determined if p < 0.01 using the paired Student’s t-Test between the two distributions.
'''


cutoff = {'pval': 0.01, 'cn': 0.5}
#cutoff = {'pval': 0.05, 'cn': 1} #test

def get_loss_allels(loh, cutoff):
    if not os.path.getsize(loh):
        return []
    df = pd.read_table(loh)
    df = df.drop_duplicates()
    df = df[(df.PVal_unique < cutoff['pval'])] #Allelic imbalance is determined if p < 0.01
    df = df[(df.HLA_type1copyNum_withBAFBin < cutoff['cn']) | (df.HLA_type2copyNum_withBAFBin < cutoff['cn'])] #copy number < 0.5
    return df.LossAllele

def main(infile, outfile):
    out_list = list(get_loss_allels(infile, cutoff))
    #sys.stdout.write(out_list + '\n')
    print(out_list)
    with open(outfile, 'w') as fw:
        for x in out_list:
            fw.write(x + '\n')

if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) != 2:
        #sys.stdout.write('usage: python {} infile(lohhla xls) outfile(result)\n'.format(sys.argv[0]))
        print('usage: python {} infile(lohhla xls) outfile(result)'.format(sys.argv[0]))
        exit(1)
    infile, outfile = args
    main(infile, outfile)
