import gzip
import sys
helper = "call ~ <fname1> <fname2> <fname3> ...\n\t\ttransfrom the gzipped fastq file to fasta format, ignore the base quality line."
if len(sys.argv) == 1:
    print(helper)
    print("exit!")
    exit()
    
for fname in sys.argv[1:]:
    if len(fname[-2:]) == 'gz':
        cnt = 1
        with gzip.open(fname) as fin, with open(fname[:fname.find('.')]+'.fa', 'w') as fout:
            for l in fin:
                if cnt % 4 == 1:
                    ans = ">"
                    for c in l[1:-1]:
                        if c.isalnum():
                            ans += c
                        else:
                            ans += '_'
                    fout.write(ans)
                    fout.write('\n')
                elif cnt % 4 == 2:
                    fout.write(l)
                cnt += 1
    else:
        print(f"{fname} not ending as gz, ignored.\nPlease verify that it's a gzipped fastq file!")
