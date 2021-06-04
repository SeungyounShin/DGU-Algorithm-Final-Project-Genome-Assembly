f = open('ref.txt')
ref = f.read()
f.close()

f = open('recon_vec.txt')
recon = f.read()
f.close()

recon = recon.split('\n')[:-1]

total_match = 0
for r in recon:
    best_matched = 0
    rlen = len(r)
    for i in range(len(ref)-rlen):
        matched = 0
        refSub = ref[i:i+rlen]
        for j in range(rlen):
            if(refSub[j]==r[j]):
                matched += 1
        if(best_matched < matched):
            best_matched = matched
    total_match += best_matched

print(total_match, len(ref))
