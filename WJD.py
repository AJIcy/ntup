#Skeleton joboption for a simple analysis job

include( "WJDNtup/configureServices.py" )

#---- Minimal job options -----

jps.AthenaCommonFlags.AccessMode = "ClassAccess"              #Choose from TreeAccess,BranchAccess,ClassAccess,AthenaAccess,POOLAccess
#jps.AthenaCommonFlags.TreeName = "MyTree"                    #when using TreeAccess, must specify the input tree name
jps.AthenaCommonFlags.HistOutputs = ["MYSTREAM:tree_test7.root"]  #register output files like this. MYSTREAM is used in the code

from PerfMonComps.PerfMonFlags import jobproperties as jp
jp.PerfMonFlags.doMonitoring = False
jp.PerfMonFlags.doFastMon = False


# ToolSvc += CfgMgr.GoodRunsListSelectionTool("GRLTool",GoodRunsListVec=["data16_13TeV.periodAllYear_DetStatus-v89-pro21-01_DQDefects-00-02-04_PHYS_StandardGRL_All_Good_25ns.xml"])
# filterSeq = CfgMgr.AthSequencer("AthFilterSeq")
# filterSeq += CfgMgr.GRLSelectorAlg(Tool=ToolSvc.GRLTool)

# svcMgr.MetaDataSvc.MetaDataTools += [ "LumiBlockMetaDataTool" ]
# theApp.CreateSvc += ['xAOD2NtupLumiSvc']



athAlgSeq += CfgMgr.WJDNtupAlg( CascadeFit = False, #doesn't work as of 29.07.2018, there is a memory leak in the tool. It is able to process about 0.75-1 mil. Ds candidates before shutting down.
							doTruth = False, #Do "Truth" plots, very slow.
							doData = True, #Do Data plots
							doNtup = True,
							Verbose = False,
							DsIsoCheck = False,
							SeedTrackCut = False,
							Iso_cut = False,
							Collinear_cut = 0.8,
							Jpsi_range = 6, #mm. range from J/psi to loop over tracks for SeedTrackCut try 20, 10, 6, and W/O or 300
							SeedpTMin = 1.5, #GeV
							SeedpTMax = 1000, #GeV
							SeedJpsidRMin = 1.5, #Pi, seed track V Jpsi track
							dRSeedCone = 1, #Pi, dR of secondary tracks in relation to the seed track, only has meaning if "SeedTrackCut" is true
							JpsiPtCut = 10, #GeV
							DsPtCut = 10, #GeV
							VertexFitter = getVKalVertexFitter(), #Fitter to test the first two candidates for a Ds
							JpsiVertexFitter = getVKalVertexFitterJ(), #J/psi fitter
							DsVertexFitter = getVKalVertexFitterDs(), #Ds fitter
							TrackSelectionTool = TrackSelectionTool)                              #adds an instance of your alg to the main alg sequence


#---- Options you could specify on command line -----
#jps.AthenaCommonFlags.EvtMax=50                          #set on command-line with: --evtMax=-1
#jps.AthenaCommonFlags.SkipEvents=0                       #set on command-line with: --skipEvents=0
#jps.AthenaCommonFlags.FilesInput = ["/eos/user/v/vtskhaya/WJD/WJD_AOD_1/user.vtskhaya.WJD2.AOD.merge2.pool.v1_EXT0/user.vtskhaya.14241422.EXT0._000001.AOD.merge.pool.root","/eos/user/v/vtskhaya/WJD/WJD_AOD_1/user.vtskhaya.WJD2.AOD.merge2.pool.v1_EXT0/user.vtskhaya.14241422.EXT0._000003.AOD.merge.pool.root","/eos/user/v/vtskhaya/WJD/WJD_AOD_1/user.vtskhaya.WJD2.AOD.merge2.pool.v1_EXT0/user.vtskhaya.14241422.EXT0._000004.AOD.merge.pool.root","/eos/user/v/vtskhaya/WJD/WJD_AOD_1/user.vtskhaya.WJD2.AOD.merge2.pool.v1_EXT0/user.vtskhaya.14241422.EXT0._000005.AOD.merge.pool.root","/eos/user/v/vtskhaya/WJD/WJD_AOD_1/user.vtskhaya.WJD2.AOD.merge2.pool.v1_EXT0/user.vtskhaya.14241422.EXT0._000006.AOD.merge.pool.root","/eos/user/v/vtskhaya/WJD/WJD_AOD_1/user.vtskhaya.WJD2.AOD.merge2.pool.v1_EXT0/user.vtskhaya.14241422.EXT0._000007.AOD.merge.pool.root","/eos/user/v/vtskhaya/WJD/WJD_AOD_1/user.vtskhaya.WJD2.AOD.merge2.pool.v1_EXT0/user.vtskhaya.14241422.EXT0._000008.AOD.merge.pool.root","/eos/user/v/vtskhaya/WJD/WJD_AOD_1/user.vtskhaya.WJD2.AOD.merge2.pool.v1_EXT0/user.vtskhaya.14241422.EXT0._000009.AOD.merge.pool.root","/eos/user/v/vtskhaya/WJD/WJD_AOD_1/user.vtskhaya.WJD2.AOD.merge2.pool.v1_EXT0/user.vtskhaya.14241422.EXT0._000010.AOD.merge.pool.root"]        #set on command-line with: --filesInput=...


include("AthAnalysisBaseComps/SuppressLogging.py")              #Optional include to suppress as much athena output as possible. Keep at bottom of joboptions so that it doesn't suppress the logging of the things you have configured above
