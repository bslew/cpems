#!/usr/bin/env python2.7

'''
Module description: a bunch of common usage python functions


Created on Mar 1, 2012
@author: blew
'''

import sys
import os
import errno
sys.path.append('./')  # FOR EXTRA MODULES LOCATED IN LOCAL DIRECTORY
from subprocess import Popen, PIPE
import subprocess

#from pylab import *
import numpy as np
from datetime import datetime, date, time
import signal
import h5py
from astropy.time import Time

import csv

#import PySide.QtCore.QSettings


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# FNUCTION SPACE

def getStdOutLine(cmd,searchKey):
    for line in os.popen(cmd).readlines():     # run find command
        if line.find(searchKey)!= -1:
            return line.strip()

    return ""

def getStdOutValue(cmd,searchKey, separator=":"):
    for line in os.popen(cmd).readlines():     # run find command
        if line.find(searchKey)!= -1:
            val=line.strip().split(separator)
            return val[-1].strip()
    return ""

def getStdOutValues(cmd,searchKeys, separator=":"):
    vals=dict()
    for line in os.popen(cmd).readlines():     # run find command
        for key in searchKeys:
            searchKey=key
            if line.find(searchKey)!= -1:
                val=line.strip().split(separator,1)
                vals[searchKey]=val[-1]
    return vals

#def isArray(var):
#    return lambda var: isinstance(var, (list, tuple))
#    if type(var)==type(list()): return True
#    return False

def isList(var):
#    return lambda var: isinstance(var, (list))
    if type(var)==type(list()): return True
    return False


def readTxtFileAsStringArray(fname,separator=" "):
    f=open(fname)
    lines=f.readlines()
    a=[]
    for l in lines:
        a.append(l.strip().split(separator))
    return a

def executeCommandStr(cmd, printOutput=False):
    if printOutput:
        os.system(cmd)
    else:
#        proc = Popen([os.environ["MSCS_PROGRAM_DIR"]+'/bin/'+cmd, ''], stdout=PIPE, shell='/bin/bash')
        proc = Popen([cmd, ''], stdout=PIPE, shell='/bin/bash')
#    output = Popen([os.environ["MSCS_PROGRAM_DIR"]+'/bin/'+cmd, ''], stdout=PIPE, shell='/bin/bash')
        print proc.stdout.read()
#    output = proc.communicate()[0]
#        print output
#    return output
#    result = os.popen(cmd).readlines()
#    if printOutput:
#        for l in result:
#            print "   ***|"+l
#    return result

def fileExists(fname):
    return os.path.isfile(fname)

def sayAndExecute(msg,cmd='',execute=0, quiet=False):
    print "----------------------------------"
    print " * "+msg
    print
    ret=0
    if cmd!='':
        print " executing command: "
        print cmd
        print
    if quiet:
        cmd+=' > /dev/null'
    if execute==1:
        cmd='/bin/bash -c "%s"' % cmd
#         ret=os.system(cmd);
#         ret = subprocess.check_call(cmd, shell=True)
#         ret = subprocess.call(cmd, shell=True)
#         return ret
#         ret = subprocess.Popen(cmd, stdout=subprocess.STDOUT, stderr=subprocess.STDOUT, shell=True)
#         ret = subprocess.Popen([cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
#         print ret.communicate()[0]
#     return ret.returncode
        try:
            cmdout = subprocess.check_call(cmd, shell=True)
            
        except subprocess.CalledProcessError as cmdexc:
            print cmdexc.returncode
            return cmdexc.returncode  
#         except KeyboardInterrupt:
#             return -1
    return ret
    
def strToFloat(v):
    if v=='':
        return 0 
    return float(v)

def strToInt(v):
    if v=='':
        return 0 
    return int(v)

def getIntList(strData):
    tmp=[ strToInt(v) for v in strData.split(',')]
    return tmp

def getFloatList(strData):
    tmp=[ strToFloat(v) for v in strData.split(',')]
    return tmp

def getCMDline():
    cmdLine=''
    for a in sys.argv:
        cmdLine=cmdLine+' '+a
    return cmdLine

def saveHowWeWereCalled():
    cmdLine=getCMDline()
    fname='.'+os.path.basename(sys.argv[0])+'.log'
    f = open(fname, "a")
    dt=datetime.now()
    f.write(dt.ctime())
    f.write("\n")
    f.write(r'%s' % cmdLine) 
    f.write("\n")
    f.close()


class AlarmException(Exception):
    pass

def alarmHandler(signum, frame):
    raise AlarmException

def waitEnter():
    print
    print
    raw_input('Press Enter to continue...')
    print
    print

# def cpeds_kbhit():
#     import select
#     dr,dw,de = select([sys.stdin], [], [], 0)
#     return dr <> []

def readUserValue(queryStr, inputType, defaultAnswer, timeoutTime=0):
    signal.signal(signal.SIGALRM, alarmHandler)
    signal.alarm(timeoutTime)
#     s = signal.signal(signal.SIGINT, signal.SIG_IGN)
#     signal.signal(signal.SIGINT, s)


    ansOut=0
    errCode=0
    if inputType=='float':
        ansOut=0.0
    if inputType=='int':
        ansOut=0
    if inputType=='string':
        ansOut=''
    print
    print 'DEFAULT ANSWER: ',defaultAnswer
    
    if timeoutTime>0:
        queryStr=queryStr+'(The DEFAULT ANSWER  will fire in %i seconds).\nYour answer' % ( timeoutTime)
#    queryStr+=' (default answer is: ',defaultAnswer
    queryStr+=':\n> '

    try:
        ans=raw_input(queryStr)
        signal.alarm(0)
    except AlarmException: 
        # timeout
#         print 'AlarmException, continuing...'
        print "using default answer:", defaultAnswer
        return defaultAnswer,0
    except EOFError: 
        # exception for script mode executions
        return defaultAnswer,0
#     except KeyboardInterrupt:
#         signal.alarm(0)
#         print 'Alarm cancelled'
#         pass
#        return ansOut,errCode
    
    if ans=='':
        ans=defaultAnswer

    try:
        if inputType=='float':
            ansOut=float(ans)
        if inputType=='int':
            ansOut=int(ans)
        if inputType=='string':
            ansOut=ans
    except ValueError:
        errCode=1
    
    return ansOut,errCode

def QsaveParameter(paramName,paramValue, quiet=False):
    from PySide.QtCore import QSettings

    if quiet==False:
        print 'saving Qparameter: %s' % paramName
        print '    parameter value: ',paramValue
    settings = QSettings("TCfA", "CPEDS")
    settings.setValue(paramName, paramValue)

def QloadParameter(paramName, quiet=False):
    from PySide.QtCore import QSettings
    settings = QSettings("TCfA", "CPEDS")
    paramValue=settings.value(paramName)
    
    if quiet==False:
        print 'loading Qparameter: %s' % paramName
        print '    parameter value: ',paramValue
    return paramValue

def saveParameter(fname,paramName,paramValue):
    f = open(fname, "a")
    f.write(r'#param %s: %s\n' % (paramName,paramValue)) 
    f.close()

def loadParameter(fname,paramName):
    params=readTxtFileAsStringArray(fname)
#    f = open(fname, "a")
#    f.write(r'#param %s: %s\n' % (paramName)) 
#    f.close()

def chdir(dirName):
    os.chdir(dirName)
    print
    print 'Changing working directory to: ',dirName
    print



def get_ij_max2Darray(a):
    Nx=len(a)
    Ny=len(a[0])
    maxval=a[0][0]
    ij=(0,0)
    for i in np.arange(Nx):
        for j in np.arange(Ny):
            if a[i][j]>maxval:
                maxval=a[i][j]
                ij=(i,j)
    return ij


def get_ij_min2Darray(a):
    Nx=len(a)
    Ny=len(a[0])
    minval=a[0][0]
    ij=(0,0)
    for i in np.arange(Nx):
        for j in np.arange(Ny):
            if a[i][j]<minval:
                minval=a[i][j]
                ij=(i,j)
    return ij


def verifyZeroExitStatusCode(code):
    if code!=0:
        print 'last command exit status code was non-zero, will exit now'
        sys.exit(code)


def saveTextToFile(fname,text):
    f=open(fname,'w')
    f.write(text)
    f.close()
    

def getHDF5key(fname,dset,key):
#    print fname
#    print dset 
#    print key
    f = h5py.File(fname, 'r+')
    dset = f['/'+dset]
    k=dset.attrs[key]
    f.close()
    return k

def getHDF5dset(fname,dset):
    f = h5py.File(fname, 'r+')
    dset = f['/'+dset]
    return dset

def h5f_node_exists(f,node):
    e = True
    try:
        f[node]
    except KeyError:
        e = False # now we know it doesn't    
    return e


def say(text,verbocity=0):
    if verbocity==0:
        print '-----------------------------------------'
        print '* ',text.upper()
        print '-----------------------------------------'


def getUDPdatagram(ip,port,N, multicast=False):
    import socket
    import struct
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM, socket.IPPROTO_UDP)
    s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
#     if multicast:
#         s.setsockopt(socket.IPPROTO_IP, socket.IP_MULTICAST_TTL, 32) 
#         s.setsockopt(socket.IPPROTO_IP, socket.IP_MULTICAST_LOOP, 1)
    
    s.bind((ip, port))
    if multicast:
#         print "subscribing to multicast"
        MCAST_GRP=ip
#         host = socket.gethostbyname(socket.gethostname())
        mreq = struct.pack("4sl", socket.inet_aton(MCAST_GRP), socket.INADDR_ANY)
        s.setsockopt(socket.IPPROTO_IP, socket.IP_ADD_MEMBERSHIP, mreq)
#         s.setsockopt(socket.SOL_IP, socket.IP_MULTICAST_IF, socket.inet_aton(host))
#         s.setsockopt(socket.SOL_IP, socket.IP_ADD_MEMBERSHIP, socket.inet_aton(MCAST_GRP) + socket.inet_aton(host))
    print "waiting on port:", port
    i=0
    allData=list()
    while 1:
        data, addr = s.recvfrom(1500)
        allData.append(data)
        i+=1
        if i==N:
            return allData
        
        

def jd2cal(jd,precision=0):
    t=Time(jd, format='jd',precision=precision)
#     t.format('%Y-%m-%d %H:%M:%S')
#     t.format='unix'
#     t.format(format)value=int(t.value)
    t.format='iso'
#     print t.jd1,t.jd2
#     print t.value
    return t.value

def u2cal(utm,precision=0):
    t=Time(utm, format='unix',precision=precision)
#     t.format('%Y-%m-%d %H:%M:%S')
#     t.format='unix'
#     t.format(format)value=int(t.value)
    t.format='iso'
#     print t.jd1,t.jd2
#     print t.value
    return t.value

def cal2u(date_time_str):
    t=Time(date_time_str, format='iso',scale='utc')
#     t.format('%Y-%m-%d %H:%M:%S')
    t.format='unix'
#     t.format(format)value=int(t.value)
#     t.format='iso'
#     print t.jd1,t.jd2
#     print t.value
    return t.value

'''
date_time_str in format %Y-%m-%s %H:%M:%S 
'''
def cal2jd(date_time_str):
    t=Time(date_time_str, format='iso',scale='utc')
#     t.format('%Y-%m-%d %H:%M:%S')
#     t.format='unix'
#     t.format(format)value=int(t.value)
    t.format='jd'
    
#     print date_time_str,t.jd1,t.jd2,t.jd
#     print t.value
    return t.jd

'''
time_interval - float [JD]
'''
def addTimeToDate(date_time_str,time_interval):
    jd=cal2jd(date_time_str)
    jd=jd+time_interval
    return jd2cal(jd)


################################################################################################################################################
'''
calculates polynomial values for arguments given by X 
a - polynomial coefficients.
returns array Y of size X with values
'''
def mkPolynomial(X,a):
    Y=list()
    for x in X:
        f=0
        for i in range(len(a)):
            f=f+a[i]*pow(x,i)
        Y.append([x,f])
    Y=np.asarray(Y)
    return Y

################################################################################################################################################
def mkDir(dirName):
    if os.path.isdir(dirName)==False:
        os.mkdir(dirName)
        
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise        

################################################################################################################################################
'''
namePattern - eg. .txt. The namePattern is used to search files in dirName
dirName - directory to search in
maxModifTimeInterval - time interval in hours. If given then the resulting list will contain only files
    such that there exist a pair for which modification time difference is not bigger than maxModifTimeInterval
    
refmtime - reference modification time object (returned by os.path.getmtime). 
    If given then only the files that yield this criteria will be on the list
    otherwise, the first files from the dirName will be taken according in order of sorted file names
'''
def getFilesList(namePattern='',dirName='.', maxModifTimeInterval=None, refmtime=None):
    l=list()
    if dirName=='':
        dirName='./'
    
    filesListAll=os.listdir(dirName)
    filesList=list()
    if refmtime==None:
        filesList=filesListAll
    else:
        for file in filesListAll:
            if (np.abs(os.path.getmtime(dirName+"/"+file)-refmtime))/3600 < maxModifTimeInterval:
                filesList.append(file)
        
    for file in sorted(filesList):
#         print file
        if namePattern in file:
            if len(l)==0:
                l.append(file)
            else:
                if maxModifTimeInterval==None:
                    l.append(file)
                else:
#                     print os.path.getmtime(file)
                    if (np.abs(os.path.getmtime(dirName+"/"+file)-os.path.getmtime(dirName+"/"+l[-1])))/3600 < maxModifTimeInterval:
                        l.append(file)

#     fl=sorted(l)
#     for file in fl:
#         if (np.abs(os.path.getmtime(file)-os.path.getmtime(l[-1])))/3600 < maxModifTimeInterval:
#             l.append(file)
        

#             print(os.path.join("/mydir", file))
    return l

################################################################################################################################################
'''
namePattern - eg. .txt. The namePattern is used to search files in dirName
dirName - directory to search in
maxModifTimeInterval - time interval in hours. If given then the resulting list will contain only files
    such that there exist a pair for which modification time difference is not bigger than maxModifTimeInterval
    
mtime_st/en - reference modification time object (returned by os.path.getmtime) that indicates initial/final time
'''
def getFilesList_mtimeRange(namePattern,dirName, mtime_st, mtime_en):
    l=list()
    if dirName=='':
        dirName='./'
    
    filesListAll=os.listdir(dirName)
    filesList=list()
    for file in filesListAll:
        if os.path.getmtime(dirName+"/"+file)>mtime_st and os.path.getmtime(dirName+"/"+file)<mtime_en:
            filesList.append(file)
        
    for file in sorted(filesList):
        if namePattern in file:
            l.append(file)

    return l

################################################################################################################################################
'''
fname - existing file name which will be used as modification time reference eg. ABC.txt. If the file is not in
    the current directory then full path should be provided.
dirName - directory to search in
maxModifTimeInterval - time interval in hours. If given then the resulting list will contain only files
    such that there exist a pair for which modification time difference is not bigger than maxModifTimeInterval
'''
def getFilesModifiedRel(fname,dirName, maxModifTimeInterval):
    l=list()
    if dirName=='':
        dirName='./'
    for file in sorted(os.listdir(dirName)):
        if (np.abs(os.path.getmtime(dirName+"/"+file)-os.path.getmtime(fname)))/3600 < maxModifTimeInterval:
            l.append(file)

    return l
    

################################################################################################################################################
def removeOutliersByMinimizingSampleVariance(data,col,thres, Nmin=1):
    from pyCPEDScommonFunctions import OutliersMinVar
    
    out=OutliersMinVar.OutliersMinVar(data=data,col=col,Nmin=Nmin)
    out.findOutliers(thres)
    
    return out.getCleanData(),out.shist
    
################################################################################################################################################
def removeOutliersByMinimizingSampleMAD(data,col,thres, Nmin=1):
    from pyCPEDScommonFunctions import OutliersMinVar
    
    out=OutliersMinVar.OutliersMinMAD(data=data,col=col,Nmin=Nmin)
    out.findOutliers(thres)
    
    return out.getCleanData(),out.shist
    
    
def mad(data, axis=None):
    return np.median(np.abs(data - np.median(data, axis)), axis)
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# MAIN PROGRAM




################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################


