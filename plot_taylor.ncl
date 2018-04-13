;*************************************************
; plot_taylor.ncl
;*************************************************
;
; Taylor plot for ETOOL package
; additional program for the ETOOL-VIS-package
; for producing Taylor plots (either temporally
;  or spatially) for monthly ETOOL output of
;  COSMO-CLM w.r.t. E-OBS for the different 
;  Prudence regions
; input parameters:
;  see user input section below
; needs: 
;  - input files as prepared by ETOOL
;  - Function 'taylor_diagram.ncl', available at:
;    https://www.ncl.ucar.edu/Applications/taylor.shtml
; call as:
;  ncl plot_taylor.ncl
;
; Susanne Brienen, DWD, 26/02/2015
;*************************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "taylor_diagram.ncl"

begin

  ;-------------------------------------------------------------
  ; USER INPUT SECTION
  ;-------------------------------------------------------------
 inpath="/scratch/snx3000/ssilje/E-Tool-out/"
  ;outpath="/scratch/snx3000/ssilje/E-Tool-out/plots/"
  year1=1981  ; start year
  year2=2000  ; end year
;  runids=(/"CON611"/)
; only needed if not called by RUN_compare.sh
; and transferred with the call of this ncl-routine:
;  var="T_2M" ; variable name
;  runids=(/"CON502_DWD","CON516","CON534","CON537","CON554"/) ; vector of experiment IDs to be compared
;  vartype="T" ; S: spatial variability, T: temporal variability
;  outpath="/pool/data/CCLM-EVAL/ANA_CON031/plots/"
;  outpath1="/pool/data/CCLM-EVAL/ANA_CON031/plots/epsplots/"
  ;-------------------------------------------------------------
  ; End USER INPUT SECTION
  ;-------------------------------------------------------------

  ; Prudence regions
  regs   = (/"BI","IP","FR","ME","SC","AL","MD","EA"/)
  nregs  = dimsizes(regs)

  ; runs
  ncases = dimsizes(runids)
  str_runids =""
  do rs = 0, ncases - 1
    str_runids = str_runids  + "_" + runids(rs)
  end do

  ; time periods
  periods = (/"DJF","MAM","JJA","SON","year"/)
  np = dimsizes(periods)

  ; different file name for TOT_PREC
  if (var .eq. "TOT_PREC") then
      nameproc="monsum"
  else
      nameproc="monmean"
  end if

  ; differences for spatial and temporal variability
  if (vartype .eq. "S") then
    ttext = "Spatial"  ; text part for plot titles
  else if (vartype .eq. "T") then
    ttext = "Temporal"
  else
    print("error in variable 'vartype'")
  end if
  end if


  ; arrays to be passed to taylor plot
  ratios = new ((/ncases, nregs, np/),float)
  corrs  = new ((/ncases, nregs, np/),float)

  ;-------------------------------------------------------------
  ; read data                 
  ;-------------------------------------------------------------

  do jj = 0, nregs-1  ; loop over regions

     do ii = 0, ncases-1

        prefix=inpath+"ANA_"+runids(ii)+"/taylor/"+runids(ii)+"_"+var+"_"+nameproc+"_"+vartype
        fname1 = prefix+"NVAR_CLM-"+obsname+"_"+year1+"-"+year2+"_seas_REG"+(jj+1)+".nc" 
        fname2 = prefix+"NVAR_CLM-"+obsname+"_"+year1+"-"+year2+"_year_REG"+(jj+1)+".nc" 
        fname3 = prefix+"COR_CLM-"+obsname+"_"+year1+"-"+year2+"_seas_REG"+(jj+1)+".nc" 
        fname4 = prefix+"COR_CLM-"+obsname+"_"+year1+"-"+year2+"_year_REG"+(jj+1)+".nc" 

        print("read file: "+fname1)
        f1 = addfile(fname1,"r") 
        ;ratios(ii,jj,0:3) = f1->$var$
        V1 = f1->$var$
        print("read file: "+fname2)
        f2 = addfile(fname2,"r") 
        V2 = f2->$var$
        ;ratios(ii,jj,4) = f2->$var$
        print("read file: "+fname3)
        f3 = addfile(fname3,"r") 
        V3 = f3->$var$
        ;corrs(ii,jj,0:3) = f3->$var$
        print("read file: "+fname4)
        f4 = addfile(fname4,"r") 
        V4 = f4->$var$
        ;corrs(ii,jj,4) = f4->$var$

        ; Get dimensions of each variable
        V1_dims =  dimsizes(dimsizes(V1))
        V2_dims =  dimsizes(dimsizes(V2))
        V3_dims =  dimsizes(dimsizes(V3))
        V4_dims =  dimsizes(dimsizes(V4))

        ; If dimension of varible is 4 reduce it to 3
        if (V1_dims .eq. 4) then
           V1_new = V1(:,0,:,:)
           delete(V1)
           V1 = V1_new
           delete(V1_new)
        end if
        if (V2_dims .eq. 4) then
           V2_new = V2(:,0,:,:)
           delete(V2)
           V2 = V2_new
           delete(V2_new)
        end if
        if (V3_dims .eq. 4) then
           V3_new = V3(:,0,:,:)
           delete(V3)
           V3 = V3_new
           delete(V3_new)
        end if
        if (V4_dims .eq. 4) then
           V4_new = V4(:,0,:,:)
           delete(V4)
           V4 = V4_new
           delete(V4_new)
        end if

        ; fill matrix of ratios and correlations
        ratios(ii,jj,0:3) = V1
        ratios(ii,jj,4) = V2
        corrs(ii,jj,0:3) = V3
        corrs(ii,jj,4) = V4

        delete([/f1,f2,f3,f4,V1,V2,V3,V4/])
   
     end do

  end do



  ;-------------------------------------------------------------
  ; plotting
  ;-------------------------------------------------------------

  markers = new(ncases,float)
  markers(:)=16
  colors0 = (/ "blue" , "red", "green", "cyan", "black" \
                         , "brown", "yellow"/) 
  colors = colors0(0:ncases-1)

  ; for spatial variability: panel and individual plots per region

    ; plot all regions on one page
    ; ----------------------------

    res= True
    plotA=new(8,graphic)

    wks = gsn_open_wks("pdf",outpath+"taylor-"+str_lower(ttext)+"_"+var+"_compare"+str_runids)

    res@Markers      = markers              ; make all solid fill
    res@Colors       = colors

    res@varLabelsYloc = 1.4                ; Move location of variable labels [default 0.45]
    res@caseLabelsFontHeightF = 0.11       ; make slight larger   [default=0.12 ]
    res@varLabelsFontHeightF  = 0.015      ; make slight smaller  [default=0.013]
    res@tiMainFontHeightF   = 0.035
    res@txFontHeightF      = 0.023

    res@stnRad        = (/ 0.5, 1.5 /)     ; additional standard radii
    res@ccRays        = (/ 0.6, 0.9 /)     ; correllation rays
    res@centerDiffRMS = True               ; RMS 'circles'

    res@taylorDraw            = False   ; do not draw picture
    res@taylorFrame           = False   ; do not advance frame

    ; one plot per region
    do ii = 0, nregs - 1

      res@tiMainString = "region " + regs(ii)
      plotA(ii) = taylor_diagram(wks,ratios(:,ii,:),corrs(:,ii,:),res)

    end do

    ; legend ressources
    lgres = True
    lgres@lgMarkerIndexes    = markers               ; make all solid fill
    lgres@lgMarkerColors     = colors
    lgres@vpWidthF           = 0.2                   ; width of legend (NDC)
    lgres@vpHeightF          = 0.035*ncases           ; height of legend (NDC)
    lgres@lgPerimOn          = False                 ; turn off perimeter
    lgres@lgLabelFontHeightF = .08                 ; label font height
    lgres@lgItemType         = "Markers"     ; draw markers only

    ; main panel
    mpresP= True
    mpresP@gsnMaximize = True
    mpresP@gsnFrame    = False                ; don't advance frame yet
    mpresP@gsnPanelYWhiteSpacePercent = 5
    mpresP@gsnPanelXWhiteSpacePercent = 5
    mpresP@txString = ttext + " variability for " + var
    mpresP@gsnPanelCenter        = False                ; don't center

    gsn_panel(wks,plotA,(/3,3/),mpresP)

    ; add legend
    gsn_legend_ndc(wks,ncases,runids,0.65,0.28,lgres)

    ; var Labels 
    txres = True
    txres@txFontHeightF = 0.012
    txres@txJust = "CenterLeft"              ; justify to the center left
    delta_y = 0.025
    ys = 0.26

    do kk = 1,np
       gsn_text_ndc(wks,kk+" - "+periods(kk-1), 0.91,ys,txres)
       ys = ys- delta_y
    end do

    frame(wks)

    ; plot regions individually as eps-files
    ;---------------------------------------

    rres = True

    rres@Markers      = markers
    rres@Colors       = colors
    rres@varLabels    = periods
    rres@caseLabels   = runids

    rres@varLabelsYloc = 1.4                ; Move location of variable labels [default 0.45]
    rres@caseLabelsFontHeightF = 0.35       ; make slight larger   [default=0.12 ]
    rres@varLabelsFontHeightF  = 0.015      ; make slight smaller  [default=0.013]

    rres@stnRad        = (/ 0.5, 1.5 /)     ; additional standard radii
    rres@ccRays        = (/ 0.6, 0.9 /)     ; correllation rays
    rres@centerDiffRMS = True               ; RMS 'circles'

    ; loop over regions
    do jj = 0, nregs-1

       wksi = gsn_open_wks("eps",outpath1+"/taylor-"+str_lower(ttext)+"_"+var+"_compare" + str_runids + "_REG-"+regs(jj))
       rres@tiMainString  = ttext + " variability for " + var + " in region " + regs(jj) ; title
       rplot  = taylor_diagram(wksi,ratios(:,jj,:),corrs(:,jj,:),rres)

    end do


end




