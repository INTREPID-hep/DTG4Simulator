#!/usr/bin/env python
import os
print ('START')

########   YOU ONLY NEED TO FILL THE AREA BELOW   #########
########   customization  area #########
#SprintName = "step1_RAW2DIGI_L1_L1P2GT.py"
OutputFileName = "l1nano" # base of the output file name, they will be saved in res directory
OutputFolder = "/eos/cms/store/user/folguera/INTREPID/DTShowers/2024_04_23/"

filtered_datasets = ['pionplus',"pionminus","proton","muon"]   # list of datasets to be processed

# LIST ALL YOUR DATASETS:
datasets = {
    'pionplus': {
        'name': 'pionplus',
        'njobs': 100,
        'eventsperjob': 1000,
        'particle': 'pi+',
        'energy': 1000,
        'ouputdir': OutputFolder,
        'outputfile': 'DTSim_muonShowers',
        'outputfiletuple': 'DTSimNtuple_muonShowers',
    },
    'pionminus': {
        'name': 'pionminus',
        'njobs': 100,
        'eventsperjob': 1000,
        'particle': 'pi-',
        'energy': 1000,
        'ouputdir': OutputFolder,
        'outputfile': 'DTSim_muonShowers',
        'outputfiletuple': 'DTSimNtuple_muonShowers',
    },
    'proton': {
        'name': 'proton',
        'njobs': 100,
        'eventsperjob': 1000,
        'particle': 'proton',
        'energy': 1000,
        'ouputdir': OutputFolder,
        'outputfile': 'DTSim_muonShowers',
        'outputfiletuple': 'DTSimNtuple_muonShowers',
    },
    'muon': {
        'name': 'muon',
        'njobs': 200,
        'eventsperjob': 1000,
        'particle': 'mu-',
        'energy': 2000,
        'ouputdir': OutputFolder,
        'outputfile': 'DTSim_muonShowers',
        'outputfiletuple': 'DTSimNtuple_muonShowers',
    },
}
queue = "longlunch" # give bsub queue -- 8nm (8 minutes), 1nh (1 hour), 8nh, 1nd (1day), 2nd, 1nw (1 week), 2nw 

########   customization end   #########
'''import socket
import datetime
import platform
import hashlib

def generate_random_int():
    # Obtener el nombre del host, la hora actual y la información del sistema
    hostname = socket.gethostname()
    current_time = str(datetime.datetime.now())
    system_info = str(platform.uname())

    # Concatenar todas las cadenas
    combined_string = hostname + current_time + system_info

    # Generar un hash SHA256 de la cadena combinada
    hash_object = hashlib.sha256(combined_string.encode())

    # Convertir el hash a un número entero
    random_int = int(hash_object.hexdigest(), 16)

    return random_int
'''
import random
def generate_random_number():
    random_number = random.randint(1, 100)
    return random_number


print ('do not worry about folder creation:')
os.system("rm -rf exec")
os.system("rm -rf /afs/cern.ch/user/f/folguera/workdir/INTREPID/batchlogs")
os.system("mkdir exec")
os.system("mkdir /afs/cern.ch/user/f/folguera/workdir/INTREPID/batchlogs")

for dataset in filtered_datasets:
    print ("dataset: "+dataset)
    print ("njobs: "+str(datasets[dataset]['njobs']))
    print ("eventsperjob: "+str(datasets[dataset]['eventsperjob']))
    print ("particle: "+datasets[dataset]['particle'])
    print ("ouputdir: "+datasets[dataset]['ouputdir'])
    print ("outputfile: "+datasets[dataset]['outputfile'])
    print()

    # First create the output folder 
    OutDir = datasets[dataset]['ouputdir']+datasets[dataset]['name']+"/"
    if not os.path.exists(OutDir):
        os.makedirs(OutDir)
    
    path = os.getcwd()
    
    #### Creating the single jobs for each dataset... 
    njobs = int(datasets[dataset]['njobs']) 
    print("Creating %d jobs" %(njobs))
    for x in range(1, int(njobs+1)):
        ##### Creates unique mac files for each job #######
        tmpfilename="exec/job_%s_%02d.mac" %(datasets[dataset]['name'],x)
        with open(tmpfilename, 'w') as fout:
            fout.write("/run/numberOfThreads 4\n")
            fout.write("/control/cout/ignoreThreadsExcept 0\n")
            fout.write("/run/initialize\n")
            fout.write("/DTSim/generator/randomizePrimary FALSE\n")
            fout.write("/DTSim/generator/sigmaMomentum 2.\n")
            fout.write("/DTSim/generator/sigmaAngle 2.\n")
            fout.write("/random/setSeeds %d %d\n" %(x*2+11,x*101))
            fout.write("/run/verbose 1\n")
            fout.write("/run/printProgress 0\n")
            fout.write("/gun/particle "+datasets[dataset]["particle"]+"\n")
            fout.write("/DTSim/generator/momentum %4.1f GeV\n" %(datasets[dataset]["energy"]))
            fout.write("/DTSim/field/value 2. tesla\n")
            fout.write("/analysis/setFileName %s/%s_%s_%dk_%02d\n" %(OutDir,datasets[dataset]['outputfile'],datasets[dataset]['name'],int(datasets[dataset]['eventsperjob']/1000),x))
            fout.write("/analysis/ntuple/setFileName 0 %s/%s_%s_%dk_%02d\n" %(OutDir,datasets[dataset]['outputfiletuple'],datasets[dataset]['name'],int(datasets[dataset]['eventsperjob']/1000),x)) 
            fout.write("/run/beamOn "+str(datasets[dataset]['eventsperjob'])+"\n")
        os.system("chmod 755 %s" %tmpfilename)


        ##### creates jobs #######
        tmpfilename="exec/job_%s_%02d.sh" %(datasets[dataset]['name'],x)
        with open(tmpfilename, 'w') as fout:
            fout.write("#!/bin/sh\n")
            fout.write("echo\n")
            fout.write("echo\n")
            fout.write("echo 'START---------------'\n")
            fout.write("echo 'WORKDIR ' ${PWD}\n")
            fout.write("cd "+str(path)+"\n")
            fout.write("source /cvmfs/geant4.cern.ch/geant4/11.2/x86_64-el9-gcc11-optdeb/bin/geant4.sh\n")
            fout.write("./exampleDTSim ./exec/job_%s_%02d.mac\n" %(datasets[dataset]['name'],x))
            fout.write("echo 'STOP---------------'\n")
            fout.write("echo\n")
            fout.write("echo\n")
        os.system("chmod 755 %s" %tmpfilename)
            
    ###### create submit.sub file ####
    with open('submit.sub', 'w') as fout:
        fout.write("executable              = $(filename)\n")
        fout.write("arguments               = $(ClusterId)$(ProcId)\n")
        fout.write("output                  = /afs/cern.ch/user/f/folguera/workdir/INTREPID/batchlogs/$(ClusterId).$(ProcId).out\n")
        fout.write("error                   = /afs/cern.ch/user/f/folguera/workdir/INTREPID/batchlogs/$(ClusterId).$(ProcId).err\n")
        fout.write("log                     = /afs/cern.ch/user/f/folguera/workdir/INTREPID/batchlogs/$(ClusterId).log\n")
        fout.write("on_exit_remove          = (ExitBySignal == False) && (ExitCode == 0)\n")
        fout.write("max_retries             = 15\n")
        fout.write("requirements            = Machine =!= LastRemoteHost\n")
        fout.write('+JobFlavour = "%s"\n' %(queue))
        fout.write("\n")
        fout.write("queue filename matching (exec/job_%s*sh)\n" %(datasets[dataset]['name']))
    
    ###### sends bjobs ######
    os.system("echo submit.sub")
    os.system("condor_submit submit.sub")
    print()
        
print()
print( "your jobs:")
os.system("condor_q")
print()
print()
