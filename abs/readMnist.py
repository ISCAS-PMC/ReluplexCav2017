
import numpy as np
import struct
import sys

def loadImageSet(filename):
 
    binfile = open(filename, 'rb')
    buffers = binfile.read()
 
    head = struct.unpack_from('>IIII', buffers, 0) 
 
    offset = struct.calcsize('>IIII') 
    imgNum = head[1]
    width = head[2]
    height = head[3]
 
    bits = imgNum * width * height 
    bitsString = '>' + str(bits) + 'B' 
 
    imgs = struct.unpack_from(bitsString, buffers, offset)
 
    binfile.close()
    imgs = np.reshape(imgs, [imgNum, width * height]) 
 
    return imgs,head
 
 
def loadLabelSet(filename):
 
    binfile = open(filename, 'rb')
    buffers = binfile.read()
 
    head = struct.unpack_from('>II', buffers, 0) 
 
    labelNum = head[1]
    offset = struct.calcsize('>II') 
 
    numString = '>' + str(labelNum) + "B" 
    labels = struct.unpack_from(numString, buffers, offset)
 
    binfile.close()
    labels = np.reshape(labels, [labelNum])
 
    return labels,head
 
 
if __name__ == "__main__":
    file1= 'train-images.idx3-ubyte'
    file2= 'train-labels.idx1-ubyte'
 
    imgs,data_head = loadImageSet(file1)
    print('data_head:',data_head)
    print(type(imgs))
    print('imgs_array:',imgs)

    if len(sys.argv) < 2:
        ord = 767
    else:
        ord = int(sys.argv[1])

    print(np.reshape(imgs[ord,:],[28,28]))

    for line in np.reshape(imgs[ord,:],[28,28]) :
        print "".join([ "".join([" "] * (4-len(str(ele)))) + str(ele) for ele in line])
        print ""
    np.savetxt(str(ord)+".txt",np.reshape(imgs[ord,:],[28,28]),fmt='%d',)
    print('----------fengexian-----------')
 
    labels,labels_head = loadLabelSet(file2)
    print('labels_head:',labels_head)
    print(type(labels))
    print(labels[1])
