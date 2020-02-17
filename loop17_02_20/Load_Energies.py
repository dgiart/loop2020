
import numpy as np
import time
# path = "D:\\Babich\\Science\\Ostapchuk\\Loop_2019\\03_2019\\functions\\"
# import sys
# sys.path.append(path)


def k():
    # print ('!!!! k - works !!!!!!')
    #Returns k-list

    k_start=time.time()
    with open('EE.txt', 'r') as f:
        ae = []
        for el in f:
            ae.append(el)
    le = []


    for el in ae:
        el = el.replace('\n', '')
        le.append(float(el.split('\t')[0]))

    _k = np.array(le)
    print (f'Time for k={time.time()-k_start}')
    return _k



def ee():
    # returns Dict:{k:EE}
    # print ('!!!! ee - works !!!!!!')
    with open('EE.txt', 'r') as f:
        ae = []
        for el in f:
            ae.append(el)


    _ee = []

    for el in ae:
        el = el.replace('\n', '')

        _ee.append(float(el.split('\t')[1]))

    _ee = np.array(_ee)
    _ee = dict(zip(k(), _ee))
    return _ee


def ek():
    # Download of EK and forming nps
    # returns Dict:{k:EK}
    # print ('!!!! ek - works !!!!!!')
    with open('EK.txt', 'r') as f:
        ak = []
        for el in f:
            ak.append(el)


    _ek = []

    for el in ak:
        el = el.replace('\n', '')

        _ek.append(float(el.split('\t')[1]))

    _ek = np.array(_ek)
    _ek = dict(zip(k(), _ek))

    return  _ek


def dee():
    # Download of dEE and forming nps
    # returns Dict:{k:dEK}
    # print ('!!!! dee - works !!!!!!')
    with open('dEE.txt', 'r') as f:
        ade = []
        for el in f:
            ade.append(el)


    _dee = []

    for el in ade:
        el = el.replace('\n', '')

        _dee.append(float(el.split('\t')[1]))

    _dee = np.array(_dee)
    _dee = dict(zip(k(), _dee))

    return  _dee

def dek():
    # Download of dEK and forming nps
    # returns Dict:{k:dEK}
    # print ('!!!! dee - works !!!!!!')
    with open('dEK.txt', 'r') as f:
        adk = []
        for el in f:
            adk.append(el)


    _dek = []

    for el in adk:
        el = el.replace('\n', '')

        _dek.append(float(el.split('\t')[1]))

    _dek = np.array(_dek)
    _dek = dict(zip(k(), _dek))
    return _dek
