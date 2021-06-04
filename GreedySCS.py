import collections

def overlap(a, b, min_length=3):
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b’s prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

def pick_maximal_overlap(reads, k):
    reada, readb = None, None
    best_olen = 0
    # Make index
    index = collections.defaultdict(set)
    for read in reads:
        for i in range(len(read) - k + 1):
            index[read[i:i+k]].add(read)
    for r in reads:
        for o in index[r[-k:]]:
            if r != o:
                olen = overlap(r, o, k)
                if olen > best_olen:
                    reada, readb = r, o
                    best_olen = olen
    return reada, readb, best_olen

def greedy_scs(reads, k):
    read_a, read_b, olen = pick_maximal_overlap(reads, k)
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a, read_b, olen = pick_maximal_overlap(reads, k)
    return ''.join(reads)

if __name__=="__main__":
    f = open("ref_10000.txt")
    ref = f.read()
    f.close()
    f = open("shortread_30_400.txt")
    shortRead = f.read().split('\n')[:-1]
    f.close()

    for k in range(2,25):
        recon = greedy_scs(shortRead, k=k)

        #print(recon[:10], len(recon))

        best_matched=0
        rlen = len(recon)
        for i in range(len(ref)-len(recon)):
            matched = 0
            refSub = ref[i:i+rlen]
            for j in range(rlen):
                if(refSub[j]==recon[j]):
                    matched += 1
            if(best_matched < matched):
                best_matched = matched

        print(k,best_matched, len(ref))
