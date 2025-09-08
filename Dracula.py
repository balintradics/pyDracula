from __future__ import print_function

import numpy as np
import sqlite3 
import csv
from sqlite3 import OperationalError

import time
import pygsl as gsl
import pygsl._numobj as numx
import pygsl.odeiv as odeiv

class Dracula:
    
    #shared vairables
    e0 = 13.606 # eV
    h = 6.626e-34 #m^{2}kg/s
    me = 9.109e-31 #kg
    mp = 1.672e-27 #kg
    kb = 1.38e-23 #J/K, or 8.617e-5 eV/K
    eV = 1.602e-19 #J
    pi = 3.141592654
    eps0 = 8.854e-12 #F/m

    #Number of all states
    #index 0 is for the free/continuum positrons
    #self.m_Nst
    #Number of active quantum states (princ. qn only), Nactst < Nst
    #self.m_Nactst; 
    
    def __init__(self,NStates=150, positron_temp=100, positron_density=1e+14, antiproton_number=100000, time=10, magnetic_field=1, usePohlrates = False, useBlackBody=False, radiationTemp = 1, useCollisionalProcesses=True, useInitialLevelPop=False, saveGSSteps = False):
        self.m_Nst = NStates
        self.m_Nactst = self.m_Nst;
        self.m_LevPop = np.zeros(self.m_Nst)
        self.m_UserLevPop = np.zeros(self.m_Nst)
        self.m_popinit = np.zeros(self.m_Nst)
        self.m_kion = np.zeros(self.m_Nst)
        self.m_ktbr = np.zeros(self.m_Nst)
        self.m_krr = np.zeros(self.m_Nst)
        
        self.m_kcol = np.zeros([self.m_Nst,self.m_Nst])
        self.m_kcol_db = np.zeros([self.m_Nst,self.m_Nst])  #new bool*[m_Nst]; // to keep track of detailed balance
        self.m_gamrad = np.zeros([self.m_Nst,self.m_Nst])
        
        self.m_Temp = positron_temp
        self.m_Density = positron_density
        self.m_NumbP = antiproton_number #initial antiproton density
        self.m_Time = time
        self.m_Bfield = magnetic_field

        #by default we don't use Pohl's scattering rates
        self.m_UsePohlRates = usePohlrates

        #by default no Black-Body radiation
        self.m_UseBlackBody = useBlackBody
        self.m_TempRad = radiationTemp

        #by default collisional processes are on
        self.m_isOverLap = useCollisionalProcesses

        #by default no initial level pop
        self.m_useLevPop = useInitialLevelPop

        # self info for each time step
        self.m_saveGSSteps = saveGSSteps
        self.m_SaveSteps = []

        # start time - arbitrary
        self.m_tStart = 0
        
        self.lam = 6.626e-34/np.sqrt(2.0*self.pi*self.me*self.kb*self.m_Temp)
        self.potfac = 8.988E9*self.eV**2
        self.k0 = self.m_Density*self.potfac**2/(self.kb*self.m_Temp*np.sqrt(self.me*self.e0*self.eV))
        self.rat = self.e0*self.eV/(self.kb*self.m_Temp)

        # Toggle e- + AntiH usecase where Ps can be produced
        self.m_psUseCase = False
        
    #properties not functions B? 
    def SetPositronTemp(self, temp):
        self.m_Temp = temp
        
    def SetBlackBodyTemp(self, temp):
        self.m_UseBlackBody = True
        self.m_TempRad = temp
    
    def SetPositronDensity(self, density):
        self.m_Density = density
    
    def SetAntiprotonN(self, n):
        self.m_NumbP = n 
    
    def SetOverlapTime(self, time):
        self.m_Time = time
    
    def SetBField(self, bfield):
        self.m_Bfield = bfield
    
    def SetUsePohlRates(self, useit):
        self.m_UsePohlRates = useit
    
    def SetPop(self, population):
        self.m_UserLevPop = population

    def SetPsUseCase(self, s = False):
        # Usecase for e- + AntiH scattering, where Ps is produced
        self.m_psUseCase = s
        
    def SetZeroTBR(self):
        self.m_ktbr = np.zeros(self.m_Nst)

    def SetZeroRR(self):
        self.m_krr = np.zeros(self.m_Nst)

    def SetStartTime(self, t):
        self.m_tStart = t
        
    def ApplyEstabPop(self) :
        for nf in range(1,self.m_Nactst):
            summ = 0 
            for ni in range(self.m_Nactst, self.m_Nst):
                summ = summ + self.m_popinit[ni] * self.m_kcol[ni][nf]
            #prev = self.m_kcol[0][nf]
            self.m_kcol[0][nf] = self.m_kcol[0][nf]+summ
            self.m_kcol[0][0] = self.m_kcol[0][0] -summ
                
        for ni  in range(1,self.m_Nactst):
            summ = 0 
            for nf in range(self.m_Nactst, self.m_Nst):
                summ = summ + self.m_kcol[ni][nf]
            #prev = self.m_kcol[ni][0]
            self.m_kcol[ni][0] = self.m_kcol[ni][0] + summ
            self.m_kcol[ni][ni] = self.m_kcol[ni][ni] - summ
        
    def UsePohlRates(self):
        self.m_UsePohlRates = True 
        print("Using Pohl et al's scattering rates for B=0 T...")
        for j in range(1,self.m_Nst):
            x = self.rat/j**2 
            self.m_kion[j]=self.k0*11.0*np.sqrt(self.rat)/((x**2.333333)+4.38*(x**1.72)+1.32*x)
            if self.m_Temp > 10 and x > 40:
                self.m_kion[j] = 0.0
            else:
                self.m_kion[j] = self.m_kion[j] * np.exp(-1*x)
                self.m_kcol[j][j] = 0
        for j in range (10, self.m_Nactst):
            self.m_ktbr[j] = self.m_kion[j]*j**2*self.m_Density*self.lam**3*np.exp(13.606*self.eV/(j**2*self.kb*self.m_Temp));
            #print("{0}\t{1}".format(j, self.m_ktbr[j]))  
        
        for ni in range(1, self.m_Nactst):
            for nf in range(1,self.m_Nst):
                if ni != nf: 
                    ei = self.rat/ni**2
                    ef = self.rat/nf**2
                    if ei > ef :
                        eg = ei
                        el = ef
                    else:
                        eg = ef
                        el = ei 
                    self.m_kcol[ni][nf] = self.k0*((ei/eg)**2.5)*ef*np.sqrt(ef)*np.exp(el-ei)*(22.0/((eg+0.9)**2.333333) + 4.5/((eg**2.5)*((eg-el)**1.33333)))
    
    def ApplyDetailedBalance(self):
        S=0.0
        for i in range(1, self.m_Nst):
            for j in range(i, self.m_Nst):
                if(self.m_kcol[i][j] > 0.0) and (i!=j) and (self.m_kcol_db[i][j] == False) and (self.m_kcol_db[j][i] == False): 
                    S = 0.5*(i**2*np.exp(13.606*self.eV/(i**2*self.kb*self.m_Temp))*self.m_kcol[i][j] + j**2*np.exp(13.606*self.eV/(j**2*self.kb*self.m_Temp))*self.m_kcol[j][i])
                    self.m_kcol[i][j] = S/(i**2*np.exp(13.606*self.eV/(i**2*self.kb*self.m_Temp)))
                    self.m_kcol_db[i][j] = True #detailed balance applied
                    self.m_kcol[j][i] = S/(j**2*np.exp(13.606*self.eV/(j**2*self.kb*self.m_Temp)))
                    self.m_kcol_db[j][i] = True # detailed balance applied

    def AddBlackBody(self) :
        print("Adding stimulated absorption/emission due to Black-Body radiation")
        Eif=0.0
        for ni in range(1,self.m_Nst):
            for nf in range(1, self.m_Nst):
                Eif = self.e0*self.eV*np.fabs((1.0/ni**2) - (1/nf**2))
                if nf<ni:
                    self.m_gamrad[ni][nf] = self.m_gamrad[ni][nf] * (1.0+ 1.0/(np.exp(Eif/(self.kb*self.m_TempRad)) - 1))
                if nf>ni:
                    self.m_gamrad[ni][nf] = self.m_gamrad[ni][nf] + self.m_gamrad[ni][nf]*(1.0/(np.exp(Eif/(self.kb*self.m_TempRad)) - 1))
                    

    def PrintLevelPop(self):
        for i in range(0,self.m_Nactst-1):
            print("{0}\t{1}".format((i+1),self.m_LevPop[i+1]))
    
    def GetThermPop(self,n):
        kbT_eV = (self.kb * self.m_Temp)/self.eV
        p = self.m_NumbP*n**2*self.m_Density*((self.h**2/(2*self.pi*self.me*self.kb*self.m_Temp))*(3.0/2.0))*np.exp(self.e0/(n**2*kbT_eV))
        return p 
    
    def SetSQLDbScat(self, dbfile, isOverlap):
        self.m_isOverLap = isOverlap
        conn = sqlite3.connect(dbfile)
        c = conn.cursor()

        # For e+ + AntiH
        if(self.m_isOverLap and self.m_psUseCase == False):
            if(self.m_Bfield == 0):
                self.UsePohlRates()
            else:
                print("Reading in scattering rates from sqlite3 file {0} ...".format(dbfile))

                for i in range(1,self.m_Nst):
                    for j in range (0, self.m_Nst):
                        ni = i
                        nf = j
                        query_base = "select RATE from XRATES where NI = {0} and NF = {1} and BFIELD = {2} and TEMP = {3}".format(ni,nf,self.m_Bfield, self.m_Temp)
                        result = c.execute(query_base)
                        returned = result.fetchone()
                        
                        if(returned != None) and (returned[0] != 0.0):
                            #print("Result: {0}".format(returned[0]))
                            if nf==0:
                                self.m_kion[i] = self.m_Density*returned[0]
                                #print("ion [{0}]: {1}".format(i, self.m_kion[i]))
                            if nf !=0:
                                self.m_kcol[i][j] = self.m_Density*returned[0]
                        else:
                             #Use Pohl rates?
                             if j == 0:
                                x = self.rat/(i**2) ;
                                self.m_kion[i] = self.k0*11.0*np.sqrt(self.rat)/(np.power(x,2.333333)+4.38*np.power(x,1.72)+1.32*x)
                                if i > 10:
                                    x =  np.log(self.m_kion[i]) - x;
                                    self.m_kion[i] = np.exp(x);
                                else:
                                    self.m_kion[i] = 0.0
                                #print("ion[{0}]: {1}".format(i,self.m_kion[i]))
                                #print("Error: non existing result for sql query:[{0}] setting ionis rate to Pohl et al rate: {1}".format(querybase,self.m_kion[i]))
                             if(j !=0) and (j!= i) and (i>10) and (j>10):
                                 ei = self.rat/i**2
                                 ef = self.rat/j**2
                                 if(ei >ef):
                                     eg = ei
                                     el = ef
                                 else:
                                     eg = ef
                                     el = ei
                                 	  #m_kcol[i][j] = k0*pow((ei/eg),2.5)*ef*sqrt(ef)*exp(el-ei)*(22.0/pow((eg+0.9),2.333333) + 4.5/(pow(eg,2.5)*pow((eg-el),1.33333)))
                                 self.m_kcol[i][j] = self.k0*np.power((ei/eg),2.5)*ef*np.sqrt(ef)*np.exp(el-ei)*(22/np.power((eg+0.9),2.333333) + 4.5/(np.power(eg,2.5)*np.power((eg-el),1.33333)))
                             if j == i : 
                                 self.m_kcol[j][j] = 0.0

                self.ApplyDetailedBalance();
                #print("Now calculate TBR rates")
                # Calculate TBR rate from ionisation rate using Saha-Boltzmann detailed balance
                TTBR=0
                for j in range(1, self.m_Nactst):
                    self.m_ktbr[j] = self.m_kion[j]*j*j*self.m_Density*self.lam**3;
#                    print(j, self.m_kion[j], self.lam, self.m_Density)
#                    print(self.lam)
                    #print("m_ktbr[{0}] = {1}".format(j, self.m_ktbr[j]))
                    if self.m_ktbr[j] == 0 :
                        self.m_ktbr[j] = 0
                    elif self.m_ktbr[j] > 0 :
                        x = np.log(self.m_ktbr[j]) + 13.606*self.eV/(j**2*self.kb*self.m_Temp);
                        self.m_ktbr[j] = np.exp(x);
                    else :
                        self.m_ktbr[j] = 0
#                    print(j, self.m_ktbr[j])
                    TTBR = TTBR + self.m_ktbr[j];
                    #     if(x < 39.0) dum1 += ktbr[j]*(1.0-exp(-kion[j]*2.e-4))/(kion[j]*2.e-4) ;
                    # print("Total tb rate {0} /s".format(TTBR))
                for j in range (self.m_Nactst-1, self.m_Nactst):
                    self.m_kion[j] = self.m_kion[self.m_Nactst-2] + (j - (self.m_Nactst-2))*(self.m_kion[self.m_Nactst-2] - self.m_kion[self.m_Nactst-3])/((self.m_Nactst-2) - (self.m_Nactst-3));
                    self.m_ktbr[j] = self.m_kion[j]*j**2*self.m_Density*self.lam**3*np.exp(13.606*self.eV/(j**2*self.kb*self.m_Temp));  
            conv = 6.12e-15 #m^3/s

            print("Reading radiative recombination rates from sqlite3 file {0}...".format(dbfile))
            # check if database exists
            query_base = "select name from sqlite_master where type = 'table' and name = 'RRATES'"
            result = c.execute(query_base)
            returned = result.fetchone()

            if returned != None :
                for i in range(1, self.m_Nst):
                    nf=i
                    query_base = "select RATE from RRATES where NF = {0} and TEMP = {1}".format(nf,self.m_Temp)
                    result = c.execute(query_base)
                    returned = result.fetchone()
                    if(returned != None) and (returned[0] != 0.0):
                        self.m_krr[nf] = self.m_Density*conv*returned[0]
                    else:
                        self.m_krr[nf] = 0.0



        ############################### 
        # For e- + AntiH, we cannot use Pohl rates for non-existing rates! Also no TBR, no RR!
        if(self.m_isOverLap and self.m_psUseCase == True):
            if(self.m_Bfield == 0):
                self.UsePohlRates()
            else:
                print("Reading in scattering rates from sqlite3 file {0} ...".format(dbfile))

                for i in range(1,self.m_Nst):
                    for j in range (0, self.m_Nst):
                        ni = i
                        nf = j
                        query_base = "select RATE from XRATES where NI = {0} and NF = {1} and BFIELD = {2} and TEMP = {3}".format(ni,nf,self.m_Bfield, self.m_Temp)
                        result = c.execute(query_base)
                        returned = result.fetchone()
                        
                        if(returned != None) and (returned[0] != 0.0):
                            #print("Result: {0}, {1}, {2}".format(ni, nf, returned[0]))
                            #print("Result: {0}".format(returned[0]))
                            if nf==0:
                                self.m_kion[i] = self.m_Density*returned[0]
                                #print("ion [{0}]: {1}".format(i, self.m_kion[i]))
                            if nf !=0:
                                self.m_kcol[i][j] = self.m_Density*returned[0]
                        else:
                            pass

            self.ApplyDetailedBalance();
                        
        conn.close()   #Close the DB 
    
    def SetSQLDbRadDec(self, dbfile, isOverlap):
        self.m_isOverLap = isOverlap            
        conn = sqlite3.connect(dbfile)
        c = conn.cursor()
        print ("Reading in decay rates from sqlite3 file {0}...this might take a while...".format(dbfile))
        for i in range(1, self.m_Nst):
            for j in range(1,self.m_Nst):
                ni = i 
                nf = j 
                query_base = "select LI,RATE from DRATES where NI = {0} and NF = {1}".format(ni,nf)
                result = c.execute(query_base)
                returned = result.fetchall()
                nrows = len(returned)
                if nrows > 0 :
                    #ncols = len(returned[0])
                    for ii in range(0, nrows,2):
                        self.m_gamrad[ni][nf]= self.m_gamrad[ni][nf]+ returned[ii][1]*(2*returned[ii][0]+1)/(ni*ni) #This needs a serious check 
                else:
                    self.m_gamrad[ni][nf] = self.m_gamrad[ni][nf] + 0 
        conn.close()

    

    # Custom function for the e- + antiH case
    def CustomAddIonizePsRates(self, dbfile):
        # enhance ionisation rates with Ps formation rates
        conn = sqlite3.connect(dbfile)
        c = conn.cursor()
        print ("Reading in Ps rates from sqlite3 file {0}...".format(dbfile))
        # check if database exists
        query_base = "select name from sqlite_master where type = 'table' and name = 'XRATES'"
        result = c.execute(query_base)
        returned = result.fetchone()                        
        if(returned != None) and (returned[0] != 0.0):
            for i in range(1,self.m_Nst):
                for j in range (0, self.m_Nst):
                    ni = i
                    nf = j
                    query_base = "select Positronium from XRATES where NI = {0} and NF = {1} and BFIELD = {2} and TEMP = {3}".format(ni,nf,self.m_Bfield, self.m_Temp)
                    result = c.execute(query_base)
                    returned = result.fetchone()
                        
                    if(returned != None) and (returned[0] != 0.0):
                        if nf==0:
                            #print("Result: {0}, {1}, {2}".format(ni, nf, returned[0]))
                            self.m_kion[i] = self.m_kion[i] + self.m_Density*returned[0]
        
    def fun(self, t, y, mu):
        f=np.zeros(self.m_Nactst+1)
        #print('func')
        ionloss = 0.0 
        iongain = 0.0 
        for i in range (1, self.m_Nactst):
            f[i] = self.m_krr[i]*y[self.m_Nactst] + self.m_ktbr[i]*y[self.m_Nactst] - self.m_kion[i]*y[i]
            ionloss = ionloss + (self.m_krr[i]*y[self.m_Nactst] + self.m_ktbr[i]*y[self.m_Nactst])
            iongain = iongain + self.m_kion[i]*y[i]
            integ = 0 
            for j in range(1, self.m_Nactst):
                if j != i:
                    integ = integ + (self.m_kcol[j][i] + self.m_gamrad[j][i])*y[j] - (self.m_kcol[i][j] + self.m_gamrad[i][j])*y[i]
            f[i] = f[i] +  integ        
        f[self.m_Nactst] = -ionloss + iongain
        return f


    def jac(self, t, y, mu):
        #print('jack')
        dfdy=np.zeros([self.m_Nactst+1,self.m_Nactst+1])
        dfdt=np.zeros(self.m_Nactst+1)
        m = dfdy 
        integ = 0.0
        integ2 = 0.0 
        for i in range(0, self.m_Nactst+1):
            dfdt[i] = 0.0 # no explicit time dependence (for now)
            integ= 0.0
            integ2 = 0.0
            for j in range(0, self.m_Nactst+1):
                if (i==0) and (j==0):
                    m[i][j] = 0.0
                else:
                    if(i != j):
                        if(i != self.m_Nactst) and ( j == self.m_Nactst):
                            m[i][j] = (self.m_krr[i]+ self.m_ktbr[i])
                        elif (i== self.m_Nactst) and (j != self.m_Nactst):
                            m[i][j] = self.m_kion[j]
                            integ2 = integ2 + self.m_krr[j] + self.m_ktbr[j]
                        else:
                            m[i][j]= self.m_kcol[j][i] + self.m_gamrad[j][i]
                            integ = integ + self.m_gamrad[i][j] + self.m_kcol[i][j]
                    if(i == j):
                        if(i == self.m_Nactst):
                            m[i][j] = -0
                            #m[i][j] = -(self.m_krr[i] + self.m_ktbr[i])

                        #else:
                            #m[i][j] = -self.m_kcol[i][j] - self.m_gamrad[i][j] - self.m_kion[i]
            if(i==self.m_Nactst):
                m[i][i]=-integ2
            if(i != self.m_Nactst):
                m[i][i] = -integ - self.m_kion[i]
            #print(type(m))
            dfdy = m
#        for i in range(0, self.m_Nactst+1):
#            for j in range(0, self.m_Nactst+1):
#                print(i,j,dfdy[i][j])
            
        return dfdy, dfdt

        
    def Compute(self):
        print("Beginning computation...")

        y = np.zeros(self.m_Nactst+1)
        #f = np.zeros(self.m_Nactst+1)
        if self.m_useLevPop == False :
            #Init the arrays
            y[self.m_Nactst] = self.m_NumbP
        else:
            for i in range(0,self.m_Nactst):
                y[i] = self.m_UserLevPop[i]
                y[self.m_Nactst] = 0


        print("Initial level population...")
        print(y)
        dimension = self.m_Nactst+1
        stepper = odeiv.step_bsimp

        step = stepper(dimension, self.fun, self.jac) 

        control = odeiv.control_y_new(step, 1e-6, 1e-6)
        evolve  = odeiv.evolve(step, control, dimension)
    
        print("# Using stepper {0} with order {1}".format(step.name(), step.order()))
        #print("# Using Control {0}".format(control.name()))
        hstart = 1e-09
        tstart = self.m_tStart
        t1 = self.m_tStart + self.m_Time*1e-06
        t = tstart
        h = hstart
        stamp = time.time()
        print("Starting integration...")
        while t < t1:
            t, h, y = evolve.apply(t, t1, h, y)
            print("Time, Stepsize, GS lev.pop.: ", t, h, y[1])
            #print("Lev. population: " , y)
            if self.m_saveGSSteps == True:
                self.m_SaveSteps.append([float(t),float(y[1])])
            
        Tot = 0
        for i in range(0, self.m_Nactst):
            self.m_LevPop[i] = y[i]
            if(i>0):
             Tot += y[i] 
        print("Total number of bound states :{0}, number of GS antiH: {1}, burned antiprotons {2}".format(Tot, y[1], self.m_NumbP - y[self.m_Nactst]))
       
    def SaveLevelPop(self, delimiter, filename):
        print("Saving result to ascii file {0}".format(filename))
        with open(filename, 'wb') as csvfile:
            writer = csv.writer(csvfile, delimiter=delimiter)
            writer.writerow(['n','pop'])
            for i in range(0,self.m_Nactst):
                writer.writerow([i+1, self.m_LevPop[i+1]])
            csvfile.close()
    
    def SaveLevelPopCSV(self, filename) :
        print("Saving result to csv file {0}".format(filename))
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            if self.m_UseBlackBody == True:
                header = "PrincQN, {0} K_BB {1} K_ {2} /m3_ {3}T_ {4} us, Thermal-Equillibrium".format(self.m_Temp, self.m_TempRad, self.m_Density, self.m_Bfield, self.m_Time)
            else:
                header = "PrincQN, {0} K_ {1} /m3_ {2}T_ {3} us, Thermal-Equillibrium".format(self.m_Temp, self.m_Density, self.m_Bfield, self.m_Time)


            header = "PrincQN, LevPop"
#            writer.writeheader(header)
            for i in range(0, self.m_Nactst-1):
                if i < 20 and self.m_Temp <= 50:
#                    writer.writerow((i+1, self.m_LevPop[i+1], "NaN"))
                    writer.writerow((i+1, self.m_LevPop[i+1]))
                elif i < 10 and self.m_Temp > 50:
#                    writer.writerow((i+1, self.m_LevPop[i+1], "NaN"))
                    writer.writerow((i+1, self.m_LevPop[i+1]))
                else:
#                    writer.writerow((i+1, self.m_LevPop[i+1], self.GetThermPop(i+1)))
                    writer.writerow((i+1, self.m_LevPop[i+1]))
            csvfile.close()
        
    def SavePop(self, i, filename, delimiter):
        print("Saving result to file {0}".format(filename))
        with open(filename, 'wb') as csvfile:
            writer = csv.writer(csvfile, delimiter=delimiter)
            for ii in range(1, self.m_Nactst):
                writer.writerow([i, self.m_LevPop[i], self.m_Temp, self.m_Density, self.m_TempRad])


    def GetLevelPop(self) :
        return self.m_LevPop

    def GetSaveSteps(self) :
        return self.m_SaveSteps


    def GetColRate_vs_Temp(self, dbfile, ni, nf, Bfield ):
        a = []
        conn = sqlite3.connect(dbfile)
        c = conn.cursor()
        query_base = "select TEMP,RATE from XRATES where NI = {0} and NF = {1} and BFIELD = {2}".format(ni,nf,Bfield)
        result = c.execute(query_base)
        returned = result.fetchall()
        if(returned != None) and (returned[0] != 0.0):
            #print("Result: {0}".format(returned))
            a = np.array(returned)

        return a


    def GetIonRate_vs_N(self, dbfile, Bfield, T ):
        a = []
        conn = sqlite3.connect(dbfile)
        c = conn.cursor()
        query_base = "select NI, RATE from XRATES where NF = 0 and BFIELD = {0} and TEMP = {1}".format(Bfield, T)
        result = c.execute(query_base)
        returned = result.fetchall()
        if(returned != None) and (returned[0] != 0.0):
            #print("Result: {0}".format(returned))
            a = np.array(returned)

        return a
    
    def GetIonRate_vs_Temp(self, dbfile, ni, Bfield ):
        a = []
        conn = sqlite3.connect(dbfile)
        c = conn.cursor()
        query_base = "select TEMP,RATE from XRATES where NI = {0} and NF = 0 and BFIELD = {1}".format(ni,Bfield)
        result = c.execute(query_base)
        returned = result.fetchall()
        if(returned != None) and (returned[0] != 0.0):
            #print("Result: {0}".format(returned))
            a = np.array(returned)

        return a

    def GetPsRate_vs_Temp(self, dbfile, ni, Bfield ):
        a = []
        conn = sqlite3.connect(dbfile)
        c = conn.cursor()
        query_base = "select TEMP,Positronium from XRATES where NI = {0} and NF = 0 and BFIELD = {1}".format(ni,Bfield)
        result = c.execute(query_base)
        returned = result.fetchall()
        if(returned != None) and (returned[0] != 0.0):
            #print("Result: {0}".format(returned))
            a = np.array(returned)

        return a
    

    def GetPsRate_vs_N(self, dbfile, Bfield, T ):
        a = []
        conn = sqlite3.connect(dbfile)
        c = conn.cursor()
        query_base = "select NI, Positronium from XRATES where NF = 0 and BFIELD = {0} and TEMP = {1}".format(Bfield, T)
        result = c.execute(query_base)
        returned = result.fetchall()
        if(returned != None) and (returned[0] != 0.0):
            #print("Result: {0}".format(returned))
            a = np.array(returned)

        return a
