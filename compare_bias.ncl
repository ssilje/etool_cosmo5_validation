;*************************************************
; compare_bias.ncl
;*************************************************
;
; additional program for the ETOOL-VIS-package
; for comparing different simulations
; plot monthly bias of COSMO-CLM w.r.t. E-OBS
;  for different Prudence regions
; input parameters:
;  see user input section below
; needs: input files as prepared by ETOOL
;  (*DIFF_CLM-EOBS*)
; call as: 
;  ncl compare_bias.ncl
;
; Susanne Brienen, DWD, 1/10/2014
;*************************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

  ;-------------------------------------------------------------
  ; USER INPUT SECTION
  ;-------------------------------------------------------------
  inpath="/scratch/snx3000/ssilje/E-Tool-out/"
;  outpath="/scratch/snx3000/ssilje/E-Tool-out/plots/"
; only needed if outpath variables are not defined in RUN_compare.sh
; and transferred with the call of this ncl-routine
;  outpath="/pool/data/CCLM-EVAL/ANA_CON031/plots/"
;  outpath1="/pool/data/CCLM-EVAL/ANA_CON031/plots/epsplots/"
  year1=1981  ; start year
  year2=2000  ; end year
; only needed if not called by RUN_compare.sh:
;  var="TMIN_2M" ; variable name
;  runids=(/"CON611"/) ; vector of experiment IDs to be compared
  ;-------------------------------------------------------------
  ; End USER INPUT SECTION
  ;-------------------------------------------------------------

  nruns = dimsizes(runids)
  str_runids =""
  do rs = 0, nruns - 1
    str_runids = str_runids  + "_" + runids(rs)
  end do

  V_avg = new((/nruns,12,8/),float)   ; average value of variable var

;!!! LOAD FILE

  if (var .eq. "TOT_PREC") then
     filepart="_DIFF_CLM-"+obsname+"_"+var+"_daysum_monsum_"+year1+"-"+year2+"_ymonmean_"
  else
     filepart="_DIFF_CLM-"+obsname+"_"+var+"_daymean_monmean_"+year1+"-"+year2+"_ymonmean_"
  end if

  do ii = 0, nruns -1  ; runs
    
    print("read data for run: "+runids(ii))

    do jj = 0, 7  ; regions

      ; read data
      filein=inpath+"ANA_"+runids(ii)+"/"+runids(ii)+filepart+"REG"+(jj+1)+".nc"
      print(" from file: "+filein)
      f1    = addfile(filein,"r")
      V1 = f1->$var$
      V1_dims =  dimsizes(dimsizes(V1))

      ; If dimension of variable is 4 reduce it to 3 
      if (V1_dims .eq. 4) then
        V1_new = V1(:,0,:,:)
        delete(V1)
        V1 = V1_new
        delete(V1_new)
      end if

      ; compute average over grid points
      do kk=0,11   ; months
         V_avg(ii,kk,jj)=avg(V1(kk,:,:))
      end do

      delete(V1)

    end do

  end do

  time=ispan(1,12,1)
  colors_list = (/"black","red","blue","darkgreen","orange","Brown","LightSkyBlue",\
                  "purple3","magenta","seagreen2"/)
  colors=colors_list(0:nruns-1)

  regions = (/"BI","IP","FR","ME","SC","AL","MD","EA"/)

;-----------------------------------------------------------------
; plot all together on one page (as pdf)
;-----------------------------------------------------------------


;!!! PLOT MODE DESIRED
  mpres= True
  plot=new(8,graphic)
;!!! OPEN PLOTTING DEVICE
  wks = gsn_open_wks("pdf",outpath+"bias_plot_"+var+"_compare"+str_runids)

;!!! DEFINE PLOT PROPERTIES 
  mpres@xyMarkLineModes =(/"Lines"/)
  mpres@trXMinF = .5
  mpres@trXMaxF = 12.5
  mpres@xyLineColors = colors  
  mpres@xyLineThicknessF = 2
  mpres@tmXBMode = "Explicit"
  mpres@tmXBValues = (/1,2,3,4,5,6,7,8,9,10,11,12/)
  mpres@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
  if (var .eq. "TOT_PREC") then
    mpres@trYMinF  = -50.
    mpres@trYMaxF  = 50.
  else if (var .eq. "CLCT") then
    mpres@trYMinF  = -0.3
    mpres@trYMaxF  = 0.3
  else
    mpres@trYMinF  = -3.0
    mpres@trYMaxF  = 3.0
  end if
  end if
  mpres@xyDashPatterns = 0
  mpres@gsnStringFontHeightF = 0.028

  mpres@gsnYRefLine           = 0.              ; reference line   

  ; plotting specifications for different variables
  if (var .eq. "T_2M" .or. var .eq. "TMIN_2M" .or. var .eq. "TMAX_2M") then
     mpres@gsnRightString = "difference [K]"   
  else if (var .eq. "PMSL") then
     mpres@gsnRightString = "difference [hPa]"
  else if (var .eq. "TOT_PREC") then
     mpres@gsnRightString = "difference [mm/month]"
  else if (var .eq. "CLCT") then
     mpres@gsnRightString = "difference [fraction]"
  else
    print("Variable "+var+" not included in script!")
  end if
  end if
  end if
  end if

  mpres@gsnDraw  = False    ;do not draw picture
  mpres@gsnFrame = False    ;do not advance frame

  ; add legend
  lgres = True
  lgres@lgLineColors=colors
  lgres@vpWidthF           = 0.2                   ; width of legend (NDC)
  lgres@vpHeightF          = 0.035*nruns           ; height of legend (NDC)
  lgres@lgPerimOn          = False                 ; turn off perimeter
  lgres@lgLabelFontHeightF = .08                 ; label font height
  lgres@lgLineThicknessF   = 2.5
  lgres@lgDashIndexes      = new(nruns,"integer",0) ; dash indexes
  lgres@lgLineLabelStrings = new(nruns,"string","")         ; no labels in lines
  ; attention: this option gives different results on different computer platforms!!!!
  ;  use only when you are sure about the result!
  ;lgres@lgItemOrder        = ispan(nruns-1,0,1)     ; reverse order of legend items
  llabels = runids

  ; plots per region
  min1=new(8,float)
  max1=new(8,float)

  do ii = 0,7
    ; compute min/max for plot
    min1(ii)=floor(min(V_avg(:,:,ii))*10)/10
    max1(ii)=ceil(max(V_avg(:,:,ii))*10)/10
;    mpres@trYMaxF = max1(ii)
;    mpres@trYMinF = min1(ii)
    mpres@gsnLeftString = regions(ii)        
    ; plot
    plot(ii) = gsn_csm_xy(wks,time,V_avg(:,:,ii),mpres)
  end do

  ; main panel
  mpresP= True
  mpresP@gsnMaximize = True
  mpresP@gsnFrame    = False                ; don't advance frame yet
  mpresP@gsnPanelYWhiteSpacePercent = 5
  mpresP@gsnPanelXWhiteSpacePercent = 5
  mpresP@txString ="BIAS CCLM - " + obsname + " ("+var+")"
  mpresP@gsnPanelCenter        = False                ; don't center
  gsn_panel(wks,plot,(/3,3/),mpresP) 

  ; legend
  gsn_legend_ndc(wks,nruns,llabels,0.71,0.28,lgres)
  frame(wks)

;-----------------------------------------------------------------
; plot regions individually as eps-files
;-----------------------------------------------------------------

  ; regions
  do ii = 0, 7
    wksi   = gsn_open_wks("eps",outpath1+"/bias_plot_"+var+"_compare" + str_runids + "_REG-"+regions(ii))
    mpres@gsnLeftString = regions(ii) 
    mpres@gsnStringFontHeightF = 0.018
    mpres@tiMainString ="BIAS CCLM-" + obsname + " ("+var+")"
    mpres@tiDeltaF = 2.4

;    mpres@trYMaxF = max1(ii)
;    mpres@trYMinF = min1(ii)

    ; plot
    ploti=gsn_csm_xy(wksi,time,V_avg(:,:,ii),mpres)

    ; legend
    gsn_legend_ndc(wksi,nruns,llabels,0.55,0.32,lgres)

    draw(ploti)
  end do

end
