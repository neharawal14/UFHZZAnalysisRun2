#!/usr/bin/python
#-----------------------------------------------
# Latest update: 2014.09.14
#-----------------------------------------------
import sys, os, pwd, subprocess
#import sys, os, pwd
import optparse, shlex, re
import time
from time import gmtime, strftime
import math
import re
def finding_mass_value(first_name) :
    pattern = r"_M(\d+)_"
    match = re.search(pattern, first_name)
    if match:
        m_124 = "M"+match.group(1)
        print(m_124)  # Output: M-124
    else:
        print("Substring not found.")
        m_124=None
    return m_124

#define function for parsing options
def parseOptions():
    global observalbesTags, modelTags, runAllSteps

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    # input options
    parser.add_option('-t', '--tag', dest='TAG', type='string',default='', help='tag to be appended to the results, default is an empty string')
    #parser.add_option('-name', '--name_input', dest='NAME_INPUT', type='string',default='', help='name of the work area for the sample year, default is an empty string')
    parser.add_option('-d', '--datasets', dest='DATASETS', type='string', default='datasets_Fall15_25ns_MiniAODv1.txt', help='txt file with datasets to run over')
    parser.add_option('-c', '--cfg', dest='CONFIGFILE', type='string', default='/scratch/osg/dsperka/Run2/HZZ4l/CMSSW_7_6_4/src/UFHZZAnalysisRun2/UFHZZ4LAna/python/templateData_76X_cfg.py', help='configuration template')
    parser.add_option('-s', '--substring', dest='SUBSTRING', type='string', default='', help='only submit datasets with this string in the name')

    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

# define function for processing the external os commands
def processCmd(cmd, quite = 0):
    #    print cmd
    #result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
#    try:
#        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT, universal_newlines=True)
#        status = 0  # If check_output succeeds, the status is 0
#    except subprocess.CalledProcessError as e:
#        output = e.output
#        status = e.returncode
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    status = process.returncode
    output = stdout if status == 0 else stderr
    print('processing command: ',cmd)
    print('\n')
    print('Output:\n   ',output)
    print('\n')
    print('status:\n ',status)
    print('\n')

    #result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    #status = result.returncode
    #output = result.stdout
    #if (status !=0):
    #    print('Error in processing command:\n   [',cmd,']')
    #    print('Output:\n   [',output,'] \n')
    #    return "ERROR!!! "+output
    #else:
    return output

def submitAnalyzer():

    # parse the arguments and options
    global opt, args
    parseOptions()

    cfgtemplate = opt.CONFIGFILE

    # save working dir
    currentDir = os.getcwd()

    tag = opt.TAG
    #name_of_input = opt.NAME_INPUT
    #outDir='resultsAna_'+tag+'_'+name_of_input
    outDir='resultsAna_'+tag
    print(" output dir:",outDir)
    if (not os.path.isdir(outDir)):
        cmd = 'mkdir '+outDir
        print(" first command here :",cmd)
        processCmd(cmd)
        print(" second command here :",cmd)
        cmd = 'mkdir '+outDir+'/cfg/'
        processCmd(cmd)

    # get the datasets
    print('[Gathering Dataset Information]')
    datasets = []
    cross_section = {}
    nfiles = {}
    nevents = {}
    datasetfiles = {}

    with open(opt.DATASETS, "r") as datasetfile:
        for line in datasetfile:

            if (line.startswith('#')): continue

            if ( not (opt.SUBSTRING=="")):
                if (not (opt.SUBSTRING in line)): continue

            dataset = line.split()[0]
            dataset = dataset.rstrip()
            dataset = dataset.lstrip()

            datasets.append(dataset)
            cross_section[dataset] = float(line.split()[1])
            
            #cmd = './das_client.py --query="file dataset='+dataset+'" --limit=10 | grep ".root"'
            #output = processCmd(cmd)
            #while ('error' in output):
            #    time.sleep(1.0);
            #    output = processCmd(cmd)
            #datasetfiles[dataset] =  output.split()
            #nfiles[dataset] = len(datasetfiles[dataset])
 
            #cmd = './das_client.py --query="dataset dataset='+dataset+' | grep dataset.nevents" --limit=0'
            #output = processCmd(cmd)
            #while ('error' in output):
            #    time.sleep(1.0);
            #    output = processCmd(cmd)
            #nevents[dataset] = output

            #print dataset,'xs:',cross_section[dataset],'nfiles:',nfiles[dataset]#,'nevents:',nevents[dataset]
            print(dataset,'xs:',cross_section[dataset])


    # submit the jobs
    print('[Submitting jobs]')
    jobCount=0

    for dataset in datasets:
        
        #continue
        filename = dataset.split('/')[1]+'_'+dataset.split('/')[2]
        first_name = dataset.split('/')[1]
        print("first name printed here", first_name)
        #outdir_title = "OUTFILE_"+str(finding_mass_value(first_name))
        outdir_title = "OUTFILE"
        cfgfile = dataset.lstrip('/')
        cfgfile = cfgfile.replace('/','_')+'.py'

        cmd = 'cp '+cfgtemplate+' '+outDir+'/cfg/'+cfgfile
        print("*** printing command: *",cmd)
        output = processCmd(cmd)

        filelist = ''
        #for f in range(0,5):
        #    if ((f+1)>nfiles[dataset]): continue
        #    filelist += '"'
        #    if (os.path.isfile('/cms/data'+datasetfiles[dataset][f])):
        #        filelist += 'file:/cms/data'
        #    filelist += datasetfiles[dataset][f]
        #    filelist += '",'
        #filelist = filelist.rstrip(',')

        #cmd = "sed -i 's~DUMMYFILELIST~"+filelist+"~g' "+outDir+'/cfg/'+cfgfile
        cmd = "sed -i 's~DUMMYFILELIST~ ~g' "+outDir+'/cfg/'+cfgfile
        output = processCmd(cmd)

        filename = dataset.split('/')[1]+'_'+dataset.split('/')[2]
        if (len(filename)>99):
          #newfilename = filename.split('-PU')[0]
          newfilename = filename.split('-106X')[0]
          filename = newfilename

        cmd  = "sed -i 's~DUMMYFILENAME~"+filename+"~g' "+outDir+'/cfg/'+cfgfile
        output = processCmd(cmd)

        cmd  = "sed -i 's~DUMMYCROSSSECTION~"+str(cross_section[dataset])+"~g' "+outDir+'/cfg/'+cfgfile
        output = processCmd(cmd)

        if (('PromptReco' in cfgfile) or ('Run2015C_25ns-05Oct' in cfgfile) or ('16Dec' in cfgfile)):
            cmd = "sed -i 's~\"PAT\"~\"RECO\"~g' "+outDir+'/cfg/'+cfgfile
            output = processCmd(cmd)

        if (('PromptReco' in cfgfile)):
            cmd = "sed -i 's~80X_dataRun2_2016SeptRepro_v4~80X_dataRun2_Prompt_v14~g' "+outDir+'/cfg/'+cfgfile
            output = processCmd(cmd)

        if (('MCRUN2_74' in cfgfile)):
            cmd = "sed -i 's~slimmedAddPileupInfo~addPileupInfo~g' "+outDir+'/cfg/'+cfgfile
            output = processCmd(cmd)

        if (('reHLT' in cfgfile)):
            cmd = "sed -i 's~\"HLT\"~\"HLT2\"~g' "+outDir+'/cfg/'+cfgfile
            output = processCmd(cmd)

        if ( ('Single' in cfgfile) and ('Run2016' in cfgfile) ):
            cmd = "sed -i 's~checkOnlySingle = cms.untracked.bool(False)~checkOnlySingle = cms.untracked.bool(True)~g' "+outDir+'/cfg/'+cfgfile
            output = processCmd(cmd)

        crabcfgfile = 'crabConfig_'+filename+'.py'             

        cmd = 'cp crabConfig_TEMPLATE.py '+outDir+'/cfg/'+crabcfgfile
        output = processCmd(cmd)

        cmd = "sed -i 's~JOBTAG~"+tag+"~g' "+outDir+'/cfg/'+crabcfgfile
        output = processCmd(cmd)

        cmd = "sed -i 's~2022_area~"+outdir_title+"~g' "+outDir+'/cfg/'+crabcfgfile
        output = processCmd(cmd)

        cmd = "sed -i 's~CFGFILE~"+outDir+"/cfg/"+cfgfile+"~g' "+outDir+'/cfg/'+crabcfgfile
        output = processCmd(cmd)

        cmd = "sed -i 's~OUTFILENAME~"+filename+"~g' "+outDir+'/cfg/'+crabcfgfile
        output = processCmd(cmd)

        cmd = "sed -i 's~DATASETNAME~"+dataset+"~g' "+outDir+'/cfg/'+crabcfgfile
        output = processCmd(cmd)

        cmd = 'crab submit -c '+outDir+'/cfg/'+crabcfgfile
        print(cmd)     

        output = processCmd(cmd)
        if ("ERROR!!!" in output):
            #print " "
            #print " "
            #print " "
            print(" Something when wrong submitting the last dataset. You should:")
            print("     1) Remove the folder in your resultsAna directory")
            print("    2) Comment out the datasets which have been submitted in the datasets txt file")
            print("     3) Rerun the SubmitCrabJobs.py with the same arguments")
            print(" ")
            print(" ")
            print(" ")
            break
        else:
            print(output)

# run the submitAnalyzer() as main() 
if __name__ == "__main__": 
    submitAnalyzer() 
