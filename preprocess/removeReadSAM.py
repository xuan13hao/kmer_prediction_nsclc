import gzip
import sys

if len(sys.argv) < 4 or len(sys.argv) > 5:
    print("Call by\n\t\t~ <fastq> <sam> <outname>")
    exit()

fqname = sys.argv[1]
samname = sys.argv[2]
outname = sys.argv[3]
mlen = 30
if len(sys.argv) == 5:
    mlen = int(sys.argv[4])

# def getM(ss):
#     total_len = 0
#     match_len = 0
#     cnt_char = 0
#     prev = 0
#     for (i, ch) in enumerate(ss):
#         if ch.isnumeric(): continue
#         total_len += int(ss[prev:i])
#         if ch == 'M':
#             cnt_char += 1
#             match_len += int(ss[prev:i])
#         prev = i+1
#     return match_len

read_map = {}
with open(samname, 'r') as fin:
    for l in fin:
        if l[0] == '@': continue
        cont = l.split('\t')
        cigar = cont[5]
        flag = cont[1]
        name = cont[0]

        if name not in read_map:
            read_map[name] = [False, False]

        # if int(cont[1]) &0x40 == 0x40:
        #     read_map[name][0] = getM(cigar) #takeRead(cigar, mlen)
        # elif int(cont[1]) &0x80 == 0x80:
        #     read_map[name][1] = getM(cigar) #takeRead(cigar, mlen)
        # else:
        #     read_map[name][0] = getM(cigar)
        #     read_map[name][1] = mlen

print(len(read_map))

fout = open(outname, 'w')
fout.close()
cnt_map = 0
cnt = 0
buffer = ""
with gzip.open(fqname, 'r') as fin:
    for l in fin:
        l = l.decode()
        buffer += l
        cnt += 1
        if cnt == 4:
            read = buffer.split()[0][1:]

            if read in read_map:
                cnt_map += 1
            else:
                fout = open(outname, 'a')
                fout.write(buffer)
                fout.close()

            # if read not in read_map:
            #     fout = open(outname, 'a')
            #     fout.write(buffer)
            #     fout.close()
            # else:
            #     if read_map[read][0] >= mlen and read_map[read][1] >= mlen:
            #         cnt_map += 1
            #     else:
            #         fout = open(outname, 'a')
            #         fout.write(buffer)
            #         fout.close()

            buffer = ""
            cnt = 0

print(cnt_map)
