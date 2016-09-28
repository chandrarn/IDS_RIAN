% THIS IS AN IDL CODE FOR REFERENCE ON BUILDING NODES
pro repair_iccd_node,shot
mdsconnect,'landau.hit'
if shot le 21126025 then begin
  mdstcl,"edit zapmain/shot="+stremo(shot)
  mdstcl,"set def .signals"
  mdstcl,"delete node iccd /noconfirm"
  mdstcl,"add node/usage=structure .iccd"
;
  mdstcl,"set def \zapmain::top.signals.iccd"
  mdstcl,"add node/usage=signal  iccd_mon"
  mdstcl,"add tag iccd_mon iccd_mon"

  mdstcl,"write"
  mdstcl,"close"
endif

;this was added to look for pulses where the the analysis tree does not exist
mdsopen,'zapanalysis',shot,status=status
if status eq 265388041 then begin
  mdstcl,"edit zapanalysis/shot="+stremo(shot)
endif else begin
  mdstcl,"edit zapanalysis/new/shot="+stremo(shot)
endelse

;stop

mdstcl,"add node/usage=structure .iccd"
mdstcl,"set def .iccd"
mdstcl,"add node/usage=structure .spectra"
mdstcl,"set def .spectra"
mdstcl,"add node/usage=text comment"
mdstcl,"add node/usage=text enter_port"
mdstcl,"add node/usage=numeric gain"
mdstcl,"add node/usage=numeric gate"
mdstcl,"add node/usage=numeric grating"
mdstcl,"add node/usage=signal  iccd_image"
mdstcl,"add node/usage=text    view"
mdstcl,"add node/usage=numeric wavelength"
mdstcl,"add node/usage=structure .bin_param"
mdstcl,"add node/usage=structure .tele_35"
mdstcl,"add node/usage=structure .tele_90"
mdstcl,"add node/usage=structure .intensity"
;stop
mdstcl,"set def \zapanalysis::top.iccd.spectra.bin_param"
mdstcl,"add node/usage=signal row_offsets"
mdstcl,"add node/usage=signal start_bin"
mdstcl,"add node/usage=signal end_bin"
;mdstcl,"add node/usage=signal scale_fac"
;mdstcl,"add node/usage=signal inst_fwhm"
;mdstcl,"add node/usage=signal inst_func"
;stop
mdstcl,"set def \zapanalysis::top.iccd.spectra.tele_35"
mdstcl,"add node/usage=numeric trans_focus"
mdstcl,"add node/usage=numeric trans_r"
mdstcl,"add node/usage=numeric iris_diam"
mdstcl,"add node/usage=text iris_hole_se"
;stop
mdstcl,"set def \zapanalysis::top.iccd.spectra.tele_90"
mdstcl,"add node/usage=numeric trans_focus"
mdstcl,"add node/usage=numeric trans_r"
mdstcl,"add node/usage=numeric iris_diam"
mdstcl,"add node/usage=text iris_hole_se"
;stop
mdstcl,"set def \zapanalysis::top.iccd.spectra.intensity"
for i=1,9 do begin
   mdstcl,"add node/usage=signal iccd_0"+stremo(i)
   mdstcl,"add tag iccd_0"+stremo(i)+" iccd_0"+stremo(i)
   mdstcl,"add node/usage=signal iccd_0"+stremo(i)+":raw_binned"
   mdstcl,"add node/usage=numeric iccd_0"+stremo(i)+":scale_fact"
   mdstcl,"add node/usage=signal iccd_0"+stremo(i)+":inst_func"
   mdstcl,"add node/usage=numeric iccd_0"+stremo(i)+":inst_fwhm"
   mdstcl,'put iccd_0'+stremo(i)+' "\ICCD_0'+stremo(i)+':RAW_BINNED * \ICCD_0'+stremo(i)+':SCALE_FACT"'
endfor
;stop
for i=10,20 do begin
   mdstcl,"add node/usage=signal iccd_"+stremo(i)
   mdstcl,"add tag iccd_"+stremo(i)+" iccd_"+stremo(i)
   mdstcl,"add node/usage=signal iccd_"+stremo(i)+":raw_binned"
   mdstcl,"add node/usage=numeric iccd_"+stremo(i)+":scale_fact"
   mdstcl,"add node/usage=signal iccd_"+stremo(i)+":inst_func"
   mdstcl,"add node/usage=numeric iccd_"+stremo(i)+":inst_fwhm"
   mdstcl,'put iccd_'+stremo(i)+' "\ICCD_'+stremo(i)+':RAW_BINNED * \ICCD_'+stremo(i)+':SCALE_FACT"'
endfor
mdstcl,"add node/usage=signal iccd_lambda"
mdstcl,"add tag iccd_lambda iccd_lambda"
mdstcl,"add tag \zapanalysis::top.iccd iccd"
;stop
;mdstcl,"clean zapanalysis"

mdstcl,"write"
mdstcl,"close"
;stop
ss,shot
;mdsput, '.zapanalysis.iccd.spectra.tele_35:iris_diam','build_with_units(31.8,"")'
;mdsput, '.zapanalysis.iccd.spectra.tele_35:iris_hole_se','build_with_units("D","")'
;mdsput, '.zapanalysis.iccd.spectra.tele_35:trans_focus','build_with_units(1.25,"")'
;stop
;mdsput, '.zapanalysis.iccd.spectra.tele_90:iris_diam','build_with_units(22.2,"")'
;mdsput, '.zapanalysis.iccd.spectra.tele_90:iris_hole_se','build_with_units("A","")'
;mdsput, '.zapanalysis.iccd.spectra.tele_90:trans_focus','build_with_units(1.771,"")'
;stop
;if shot lt 001031000 then begin
;  mdsput, '.zapanalysis.iccd.spectra.tele_35:trans_r','build_with_units(1.90,"")'
;  mdsput, '.zapanalysis.iccd.spectra.tele_90:trans_r','build_with_units(1.94,"")'
;endif else begin
;  mdsput, '.zapanalysis.iccd.spectra.tele_35:trans_r','build_with_units(1.64,"")'
;  mdsput, '.zapanalysis.iccd.spectra.tele_90:trans_r','build_with_units(1.68,"")'
;endelse
;stop
;for i=1,9 do $
;   mdsput, "\iccd_0"+stremo(i)+" \ICCD_0"+stremo(i)+":RAW_BINNED * \ICCD_0"+stremo(i)+":SCALE_FACT"


mdsclose
;stop
end