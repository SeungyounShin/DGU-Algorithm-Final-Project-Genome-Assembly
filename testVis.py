import numpy as np
import matplotlib.pyplot as plt

f = open('ref.txt')
ref = f.read()
f.close()

f = open('recon_vec.txt')
recon = f.read()
f.close()

f = open('shortreads.txt')
read = f.read()
f.close()

recon = recon.split('\n')[:-1]
read = read.split('\n')[:-1]

M = len(read)
L = len(read[0])

#mat = np.zeros((1+M+len(recon)+2,len(ref)))
#mat[0,:] = 1
#mat[1,:] = 2
#mat[1+M+1,:] = 2

matchedstr = [0 for i in range(len(recon))]
matchedInd = [0 for i in range(len(recon))]
total_match = 0
for idx,r in enumerate(recon):
    best_matched = 0
    rlen = len(r)
    for i in range(len(ref)-rlen):
        matched = 0
        refSub = ref[i:i+rlen]
        for j in range(rlen):
            if(refSub[j]==r[j]):
                matched += 1
        if(best_matched <= matched):
            best_matched = matched
            matchedInd[idx] = i
            matchedstr[idx] = r
    total_match += best_matched

print(total_match, len(ref))

#print(ref)
#print("="*len(ref))
#for i,r in enumerate(read):
#    idx = ref.find(r)
#    mat[2+i,idx:idx+len(r)] =1
#    print(" "*idx,r)
#print("="*len(ref))

"""
for i in range(len(matchedInd)):
    print(" "*matchedInd[i],matchedstr[i])
    refsub = ref[matchedInd[i]:matchedInd[i]+len(matchedstr[i])]
    for j in range(len(refsub)):
        if(refsub[j] == matchedstr[i][j]):
            mat[3+M+i,matchedInd[i]+j] =1
        else:
            mat[3+M+i,matchedInd[i]+j] = 5

plt.imshow(mat,aspect='auto')
plt.show()
"""
