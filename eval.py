f = open('ref.txt')
ref = f.read()
f.close()

f = open('recon.txt')
recon = f.read()
f.close()

print("+ ref   : ",ref[:10],'...')
print("+ recon : ",recon[:10],'...')

best_acc = 0
ind = 0

for i in range(len(recon)):
    refSub = ref[i+1:]
    if(len(refSub)<1):
        break

    matched = 0
    l = min(len(refSub),len(recon))
    for j in range(l):
        if(refSub[j] == recon[j]):
            matched += 1

    acc = matched/len(ref)

    if(acc > best_acc):
        best_acc = acc
        ind = i

print(ref)
print(" "*ind,recon)

print(" + best acc : ", best_acc)
