;cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
;				CJONES & LEILA
;
;	This program computes the lisam index
;=============================================================================
;path2='/home/scratch104/leila/LISAM_results/Australia/'

;path2='/home/clivac/LISAM_results/Americas/'
; path2='/home/clivac/LISAM_results/Asia/'
;path2='/home/clivac/LISAM_results/Asia/india/CFSR_norain/2m10m/EOF_regular/'

path2='/home/clivac/LISAM_results/Asia/india/CFSR_rain/2m10m/EOF_regular/1979_2013/'
path2='/home/clivac/LISAM_results/Asia/india/CFSR_rain/2m10m/EOF_regular/correct_cfsr/'
; choose your variables here 
xpar      = ['!5Prec','!5u10m','!5v10m', '!5q2m','!5T2m']

;---------------------------
.run /home/leila/idl_programs/subroutines/Forest_Cannon/leap_year_map.pro

;======== input of domain (lat lon for the analysis ====


; INDIA 2:
;lat0 =  5.  ; used before  default
;latf =  45.  ; used before default
;lon0 =  60. & if(lon0 lt 0) then  lon0 = lon0 + 360.0 ; used before
;lonf =  90. & if(lonf lt 0) then  lonf = lonf + 360.0 ; used before

;lon0 =  60. & if(lon0 lt 0) then  lon0 = lon0 + 360.0
;lonf =  100. & if(lonf lt 0) then  lonf = lonf + 360.0

lat0 =  5.  
latf =  35. 
lon0 =  65. & if(lon0 lt 0) then  lon0 = lon0 + 360.0 
lonf =  85. & if(lonf lt 0) then  lonf = lonf + 360.0 

;

  latm = (lat0+latf) / 2. 
  lonm = (lon0+lonf) / 2.


;==============================================
; restoring reanalysis data  
;path0='/home/sbarc/reanalysis/ncep.pentads/'

; computing new indices:

; computing new indices:

yr0=1979
yrf=2013
nyrs=(yrf-yr0)+1


fin='/home/sbarc/arc/reanalysis/cfsr/'

;---- u10m ------------------------------------
;restore,path0+'ncep.u10m.1948.2008.IDL'
;restore,fin+'cfsr.u10m.10dg.1979.2010.IDL'
restore,fin+'cfsr.u10m.10dg.1979.2013.IDL'

  id1     = where(rlon ge lon0 and rlon le lonf,mlon)     ;
  id2     = where(rlat ge lat0 and rlat le latf,mlat)     ;
  data   = cfsr(id1,id2,*)

leap_year,data,cal_year,cal_mon,cal_day,mtot,u10m
delvar,data
;


;---- v10m
;restore,fin+'cfsr.v10m.10dg.1979.2010.IDL'
restore,fin+'cfsr.v10m.10dg.1979.2013.IDL'
data   = cfsr(id1,id2,*)
leap_year,data,cal_year,cal_mon,cal_day,mtot,v10m
delvar,data
;
;----- temp
;restore,fin+'cfsr.t2m.10dg.1979.2010.IDL'
restore,fin+'cfsr.t2m.10dg.1979.2013.IDL'

  data   = cfsr(id1,id2,*)
 leap_year,data,cal_year,cal_mon,cal_day,mtot,ta2m
 delvar,data

;------ qa2m
;restore,fin+'cfsr.q2m.10dg.1979.2010.IDL'
restore,fin+'cfsr.q2m.10dg.1979.2013.IDL'
  data   = cfsr(id1,id2,*)
leap_year,data,cal_year,cal_mon,cal_day,mtot,qa2m
delvar,data
;----- prate
;restore,fin+'cfsr.prate.10dg.1979.2010.IDL'
restore,fin+'cfsr.prate.10dg.1979.2013.IDL'
  data   = cfsr(id1,id2,*)
leap_year,data,cal_year,cal_mon,cal_day,mtot,prec
delvar,data


;------------------
  rlat = rlat(id2)
  rlon = rlon(id1)
  mlat = n_elements(rlat)
  mlon = n_elements(rlon)
;=================
;
; save eofs; save eofs
ff1=strcompress(string('lisam.trdEOF.prec.d.79.2013.CFSR.',abs(fix(lat0)),'N.',abs(fix(latf)),'N.lon',fix(lon0),'.',fix(lonf),'.IDL'),/remove_all)
; figure with eofs
fileps2=path2+'LIAM_ASIA_PC1_PC2_79.2013.CFSR_prec.'+ strcompress(string(abs(fix(lat0)),'N.',abs(fix(latf)),'N.lon',fix(lon0),'.',fix(lonf), '.eps'),/remove_all)
; amplitudes e durations - files
filetxt=path2 + 'PC1.asia.ampli.dur.prec.79.2013.'+ strcompress(string(abs(fix(lat0)),'N.',abs(fix(latf)),'N.lon',fix(lon0),'.',fix(lonf), '.'),/remove_all)

filetxt2=path2 + 'PC2.india.ampli.dur.prec.79.2013.'+ strcompress(string(abs(fix(lat0)),'N.',abs(fix(latf)),'N.lon',fix(lon0),'.',fix(lonf), '.'),/remove_all)
;indices lisam1 e lisam2
fffindx=path2+'liasia1_and_2_indx_CFSR_79.2013.prec.'+ strcompress(string(abs(fix(lat0)),'N.',abs(fix(latf)),'N.lon',fix(lon0),'.',fix(lonf), '.txt'),/remove_all)
; file with lag-correlation
fsaic=strcompress(string('lisasia_lagcor.prec.cfsr.79.2013.',abs(fix(lat0)),'N.',fix(latf),'.N.lon.',fix(lon0),'.',fix(lonf),'.eps'),/remove_all)
ftrmm=path2+'Correl_liam1_liam2_TRMM_u10.v10.prec.q2.'+ strcompress(string(abs(fix(lat0)),'N.',abs(fix(latf)),'N.lon',fix(lon0),'.',fix(lonf), '.eps'),/remove_all)

;
;==============================================
;=============================================================================
;    spatial smooting
;    NOTE:  some grid points in the seleted domain (South America)
;    have 850 mb fields that are below the surface because of the
;    Andes. so applying the spatial filter is not a bad idea
;
.run /home/cjones/mjo/idl.for.subs/sub.spatial.filter.pro
;

.run
for m = 0, mtot - 1 do begin
;print,m
    map         = prec(*,*,m)
    spatfilt, map, mlon,mlat, tmp
    prec(*,*,m) = tmp
    map         = u10m(*,*,m)
    spatfilt, map, mlon,mlat, tmp
    u10m(*,*,m) = tmp
    map         = v10m(*,*,m)
    spatfilt, map, mlon,mlat, tmp
    v10m(*,*,m) = tmp
    map         = qa2m(*,*,m)
    spatfilt, map, mlon,mlat, tmp
    qa2m(*,*,m) = tmp
    map         = ta2m(*,*,m)
    spatfilt, map, mlon,mlat, tmp
    ta2m(*,*,m) = tmp
endfor
end

moyenvert,ta2m,mean,std


loadct,5
;erase
pp= [0.1,0.1,0.4,0.8]
;tvim,mean,/scale,pos = pp, $
;   stitle='    ',lcharsize=1.1,barwidth=0.75
;   oplot,[j1,j1],[i1,i2]  & oplot,[j2,j2],[i1,i2]
;   oplot,[j1,j2],[i1,i1]  & oplot,[j1,j2],[i2,i2]
; map_set, 0,lonm,/cont,/cyl,limit=[lat0,lon0,latf,lonf],/noerase,$
;	/noborder,pos = pp,mlinethick=1.9,color=100



;=============================================================================
;	Weight grid points by sqrt(cos(lat)) before computing EOF
;

xtmp     = fltarr(mlon,mlat)                               ; array with 
  for j = 0, mlon - 1 do xtmp(j,*) = sqrt(cosd(rlat))        ; sqrt(cos lat))
  xata     = fltarr(mlon,mlat,mtot)                          ; 3D array with 
  for m = 0, MTOT - 1 do xata(*,*,m) = xtmp                  ; sqrt(cos lat))
;
.run

for m = 0, MTOT - 1 do begin
    PREC(*,*,m)   = PREC(*,*,m)  * xata
    U10M(*,*,m)   = U10M(*,*,m)  * xata
    V10M(*,*,m)   = V10M(*,*,m)  * xata
    QA2M(*,*,m)   = QA2M(*,*,m)  * xata
    TA2M(*,*,m)   = TA2M(*,*,m)  * xata
endfor
end

;=============================================================================
;       remove long term mean
;
moyenvert,prec, mean_prec, sdv
moyenvert,u10m, mean_u10m, sdv
moyenvert,v10m, mean_v10m, sdv
moyenvert,qa2m, mean_qa2m, sdv
moyenvert,ta2m, mean_ta2m, sdv
;
.run
for m = 0, mtot - 1 do begin
   prec(*,*,m) = prec(*,*,m) - mean_prec
   u10m(*,*,m) = u10m(*,*,m) - mean_u10m
   v10m(*,*,m) = v10m(*,*,m) - mean_v10m
   qa2m(*,*,m) = qa2m(*,*,m) - mean_qa2m
   ta2m(*,*,m) = ta2m(*,*,m) - mean_ta2m
endfor
end
;=============================================================================
;..........................................................................
;
;       reform data array 
; ----> this is done in space dimension because space is smaller than
;       time dimension when taking data from 1948-2007
;
;..........................................................................
;=============================================================================
.run
    nfi = 5 ; nfi = number of fields
 ;  nfi = 4
    dim = nfi * MLON * MLAT	
    AUX = fltarr(dim,mtot)
      l = 0
for i = 0, MLAT - 1   do begin
  for j = 0, MLON - 1 do begin
    AUX(l,*) = reform(PREC(j,i,*))
    l = l + 1
  endfor
endfor
for i = 0, MLAT - 1   do begin
  for j = 0, MLON - 1 do begin
    AUX(l,*) = reform(U10M(j,i,*))
    l = l + 1
  endfor
endfor
for i = 0, MLAT - 1   do begin
  for j = 0, MLON - 1 do begin
    AUX(l,*) = reform(V10M(j,i,*))
    l = l + 1
  endfor
endfor
for i = 0, MLAT - 1   do begin
  for j = 0, MLON - 1 do begin
    AUX(l,*) = reform(QA2M(j,i,*))
    l = l + 1
  endfor
endfor
for i = 0, MLAT - 1   do begin
  for j = 0, MLON - 1 do begin
    AUX(l,*) = reform(TA2M(j,i,*))
    l = l + 1
  endfor
endfor
end
;=============================================================================
;        compute covariance matrix 
;
.run
  covar = aux # transpose(aux) 	; <----- covariance matrix
  corre = fltarr(dim,dim)
for i = 0, dim-1 do begin
 for j = 0, dim - 1 do begin
    corre(j,i) = covar(j,i) / sqrt(covar(j,j) * covar(i,i))
 endfor
endfor
    covar = covar / float(mtot) ; <----- covariance matrix
end
   tvim,corre,/scale
;==========================================================================
;..........................................................................
;
;
;           compute eigenvalues and eigenvectors
;
;
;..........................................................................
;==========================================================================
;
 eval = EIGENQL(CORRE,/absolute,/double,eigenvectors=evec)

 totvar = total(eval) &  perc = eval/totvar * 100.  &  plot, perc(0:20)
;==========================================================================
;          transform EOFs to arrays
; 

ff=path2 + 'eval_correl_evec_perc.IDL'
save,file=ff,corre,eval,evec

.run
       modes = 10
    EOF_PREC = fltarr(MLON,MLAT,MODES)
    EOF_U10M = fltarr(MLON,MLAT,MODES)
    EOF_V10M = fltarr(MLON,MLAT,MODES)
    EOF_QA2M = fltarr(MLON,MLAT,MODES)
    EOF_TA2M = fltarr(MLON,MLAT,MODES)
        spat = dim / nfi
for md = 0, modes - 1 do begin
;                tmp = evec(0:spat-1,md)
;   EOF_U10M(*,*,md) = reform(tmp,mlon,mlat)
;                tmp = evec(spat:spat+spat-1,md)
;   EOF_V10M(*,*,md) = reform(tmp,mlon,mlat)
;                tmp = evec(2*spat:2*spat+spat-1,md)
;   EOF_QA2M(*,*,md) = reform(tmp,mlon,mlat)
;                tmp = evec(3*spat:3*spat+spat-1,md)
;   EOF_TA2M(*,*,md) = reform(tmp,mlon,mlat)
                tmp = evec(0:spat-1,md)
   EOF_PREC(*,*,md) = reform(tmp,mlon,mlat)
                tmp = evec(spat:spat+spat-1,md)
   EOF_U10M(*,*,md) = reform(tmp,mlon,mlat)
                tmp = evec(2*spat:2*spat+spat-1,md)
   EOF_V10M(*,*,md) = reform(tmp,mlon,mlat)
                tmp = evec(3*spat:3*spat+spat-1,md)
   EOF_QA2M(*,*,md) = reform(tmp,mlon,mlat)
                tmp = evec(4*spat:4*spat+spat-1,md)
   EOF_TA2M(*,*,md) = reform(tmp,mlon,mlat)
endfor
end
;
!p.multi=[0,2,3]   & erase
  md = 1
!p.multi(0) = 5    & tvim,EOF_PREC(*,*,md),/scale
!p.multi(0) = 4    & tvim,EOF_U10M(*,*,md),/scale
!p.multi(0) = 3    & tvim,EOF_V10M(*,*,md),/scale
!p.multi(0) = 2   & tvim,EOF_QA2M(*,*,md),/scale
!p.multi(0) = 1    & tvim,EOF_TA2M(*,*,md),/scale

;=============================================================================
;             Compute PC's
;
.run
   PCS = fltarr(mtot,modes)
for m = 0, MTOT - 1 do begin
  for md = 0, modes - 1 do begin
    pcs(m,md) = total( AUX(*,m) * EVEC(*,md)) 
  endfor
endfor
end
;=============================================================================
;    adjust signs so that positive amplitudes of the time coeffs
;    correspond to positive precipitation anomalies in South America
;
; activate only if necessary
    md = 1
    EOF_PREC(*,*,md) = -EOF_PREC(*,*,md)
    EOF_U10M(*,*,md) = -EOF_U10M(*,*,md)
    EOF_V10M(*,*,md) = -EOF_V10M(*,*,md)
    EOF_QA2M(*,*,md) = -EOF_QA2M(*,*,md)
    EOF_TA2M(*,*,md) = -EOF_TA2M(*,*,md)
          PCS(*,md) = -PCS(*,md)
;=============================================================================
;        grab first two PCs
;
!p.multi=[0,1,1]  &  erase 

  liam1 =    PCS(*,0)
  liam2 =    PCS(*,1) ; attention because this is a particular case
  liam3 =    PCS(*,2)
  perc_1 =   perc(0)
  perc_2 =   perc(1)
  perc_3 =   perc(2)


;=============================================================================
;delvar,aux,corre,coar,dim,i,j,m,md,nfi,perc,res,spat,totvar,l,pcs,modes,evec
;delvar,sdv,mean_prec,mean_u10m,mean_v10m,mean_qa2m,mean_ta2m,covar,perc,eval
;=============================================================================
;       compute correlation LISAM x data for plotting
;

.run
;   nfi = 5
   cor = fltarr(mlon,mlat,3,nfi)
for  i = 0, mlat - 1 do begin
 for j = 0, mlon - 1 do begin
  tsr = reform(prec(j,i,*))
; tsr = reform(u10m(j,i,*))
   cor(j,i,0,0) = correlate(tsr,liam1) & cor(j,i,1,0) = correlate(tsr,liam2)
   cor(j,i,2,0) = correlate(tsr,liam3) 
 
; tsr = reform(v10m(j,i,*))
   tsr = reform(u10m(j,i,*))
   cor(j,i,0,1) = correlate(tsr,liam1) & cor(j,i,1,1) = correlate(tsr,liam2)
   cor(j,i,2,1) = correlate(tsr,liam3) 
 

  ; tsr = reform(qa2m(j,i,*))
   tsr = reform(v10m(j,i,*))
   cor(j,i,0,2) = correlate(tsr,liam1) & cor(j,i,1,2) = correlate(tsr,liam2)
   cor(j,i,2,2) = correlate(tsr,liam3) 
 
 ;tsr = reform(ta2m(j,i,*))
   tsr = reform(qa2m(j,i,*))
   cor(j,i,0,3) = correlate(tsr,liam1) & cor(j,i,1,3) = correlate(tsr,liam2)
   cor(j,i,2,3) = correlate(tsr,liam3) 
 
tsr = reform(ta2m(j,i,*))
   cor(j,i,0,4) = correlate(tsr,liam1) & cor(j,i,1,4) = correlate(tsr,liam2)
   cor(j,i,2,4) = correlate(tsr,liam3) 

 
endfor
endfor
end
;

erase
;tvim,cor(*,*,1,0),/scale
;-------------------------------------------------------------
;         check
; check PC1:
window,0,xsize=800,ysize=800
;                               
!p.multi=[0,2,3] & erase & k = 0
!p.multi(0)=5 & tvim,cor(*,*,k,0),/scale,range=[-1,1.,0.2],title='!5prec'
!p.multi(0)=4 & tvim,cor(*,*,k,1),/scale,range=[-1,1.,0.2],title='!5u'
!p.multi(0)=3 & tvim,cor(*,*,k,2),/scale,range=[-1,1.,0.2],title='!5v'
!p.multi(0)=2 & tvim,cor(*,*,k,3),/scale,range=[-1,1.,0.2],title='!5q2m'
!p.multi(0)=1 & tvim,cor(*,*,k,4),/scale,range=[-1,1.,0.2],title='!5T2m'
;
window,1
!p.multi=[0,1,1] & erase & k = 0
plot,pcs(*,k),xstyle=1

; check PC2:
window,2,xsize=800,ysize=800
;                               
!p.multi=[0,2,3] & erase & k = 1
!p.multi(0)=5 & tvim,cor(*,*,k,0),/scale,range=[-1,1.,0.2],title='!5prec'
!p.multi(0)=4 & tvim,cor(*,*,k,1),/scale,range=[-1,1.,0.2],title='!5u'
!p.multi(0)=3 & tvim,cor(*,*,k,2),/scale,range=[-1,1.,0.2],title='!5v'
!p.multi(0)=2 & tvim,cor(*,*,k,3),/scale,range=[-1,1.,0.2],title='!5q2m'
!p.multi(0)=1 & tvim,cor(*,*,k,4),/scale,range=[-1,1.,0.2],title='!5T2m'
;

window,3
!p.multi=[0,1,1] & erase & k = 1
plot,pcs(*,k),xstyle=1





;--------------------------------------------------------------------

fff=path2+ff1

;save,file=fff,evec,eval,perc,covar,modes,pcs,perc,nfi,lat1,latf,lon1,lonf,dim,corre,cal_mon1,cal_mon2,cal_year,cal_day1,cal_day2,eof_prec,eof_u10m,eof_v10m,eof_qa2m,eof_ta2m,liam1,liam2,perc_1,perc_2,liam3,cor,mtot,nyrs

mlat=n_elements(rlat)
mlon=n_elements(rlon)
save,file=fff,evec,eval,perc,covar,modes,pcs,perc,nfi,lat0,latf,lon0,lonf,dim,corre,cal_mon,cal_year,cal_day,eof_u10m,eof_v10m,eof_qa2m,eof_ta2m,liam1,liam2,perc_1,perc_2,liam3,cor,mtot,nyrs,rlat,rlon,mlat,mlon,xpar

;=========== computing lag correlation ===================
mlag=501
lagv=indgen(mlag)-(mlag-1)/2
coisa=fltarr(mlag)

lagcorr=c_correlate(liam1,liam2,lagv)
plot,lagv,lagcorr
oplot,lagv,coisa,line=1
res=lagcorr

id=where(abs(res) eq max(abs(res)),nid)
print,lagv(id),lagcorr(id)
; max correlation = -0.867, lag 22 pentads
id=where((res) eq max((res)),nid)
print,lagv(id),lagcorr(id) ; LISAM1 leads LISAM2 about 52 days
; max correlation = -0.867, lag 22 pentads

idmi=where(abs(res) eq min(abs(res)),nidmi)
print,lagv(idmi),res(idmi)


nada=correlate(liam1,liam2)

;----------------------------
;plotting lag-correlation and lisam-Asia and Lisam-India
!p.multi=[0,0,0]  &  erase 

fileps1=path2+fsaic
toggle,/landscape,/color,file=fileps1  


;   low  = -170     &  hig = 170  & dlt = low-30
   low  = -300     &  hig = 300  & dlt = low-40
 yr1 = min(cal_year)
  ;ref = [string(indgen(nyrs-8)+yr1 - 1900),'0'+string(indgen(8))]
;ref = [string(indgen(nyrs)+yr1 - 1900),'0'+string(indgen(8))]
ref = [string(indgen(nyrs-11)+yr1 - 1900),'0' + string(indgen(10)),'10']

  ref = strcompress(ref,/remove_all) &  print,ref
 
 vec = where(cal_day eq 1 and cal_mon eq 1)

pos = [0.05,0.1,0.95,0.75]
   plot,liam1,xrange=[0,mtot-1],xstyle=1,yrange=[low,hig],ystyle=1,yticklen=0.005, $
    xcharsize=0.0001,xticklen=0.0001, pos = pos,title=tit,charsize=1.,$
    ytitle=' ',thick=2
   oplot,[0,mtot-1],[0,0],line=1
oplot,liam2,color=150,thick=2

   for m = 0,nyrs - 1,2 do xyouts,vec(m),dlt,ref(m),alignment=0.5,charsize=1.1
   for m = 0,nyrs - 1,2 do oplot,[vec(m),vec(m)],[low,low+4]
;====
erase
;pos = [0.15,0.68,0.64,0.95]

plot,lagv,res,xstyle=1,yrange=[-1.0,1.0],ystyle=1, $
    pos = pos,charsize=1.5,charthick=2,thick=4,$
    ytitle='Correlation', xtitle='Lag (pentads)' 
nada=res
nada(*)=0
oplot,lagv,nada,line=3,thick=2

toggle
;===============================================================

;     significance level of correlations
;
ppp=0
  tsr = reform(pcs(*,ppp))
;   tsr = -1*reform(pcs(*,1))
  lag = 1
  ro1 = a_correlate(tsr,lag) & ro1 = ro1(0)
  nef = float(mtot) * (1. - ro1) / (1. + ro1)
  print,nef 
 ; tcut= 1.98
  tcut=1.96 ; for india
  sig = sqrt( (tcut^2) / (tcut^2 + nef)) & print,'pc,nef,sig',ppp+1,nef,sig
;
;
ppp=1
  tsr = reform(pcs(*,ppp))
;   tsr = -1*reform(pcs(*,1))
  lag = 1
  ro1 = a_correlate(tsr,lag) & ro1 = ro1(0)
  nef = float(mtot) * (1. - ro1) / (1. + ro1)
  print,nef 
 ; tcut= 1.98
  tcut=1.96 ; for india
  sig = sqrt( (tcut^2) / (tcut^2 + nef)) & print,'pc,nef,sig',ppp+1,nef,sig

;===========================================================
;;========= plotting using the same plot as for iteractive eof
;   plot EOFs expressed as correlations
; 
;

lat1=lat0 & lat2=latf
lon1=lon0 & lon2=lonf

toggle,/portrait,/color,file= fileps2   & loadct,0
;xpos      = fltarr(4,4)
;xpos(*,0) = [0.1,0.5,0.42,0.82] & xpos(*,1) = [0.5,0.5,0.82,0.82]
;xpos(*,2) = [0.1,0.1,0.42,0.42] & xpos(*,3) = [0.5,0.1,0.82,0.42]
;
;dx=0.22
;dy=0.40
;sp=0.07
;; left corner (bottom)
;xp1=0.05 & yp1=0.05
;;xp1=0.15 & yp1=0.1
;xp12=xp1+dx & yp12=yp1+dy
;;-------
;; left  corner (top)
;xp2=xp1     & yp2=yp12+sp
;xp22=xp2+dx & yp22=yp2+dy
;;---------
;;right corner (bottom)
;xp3=xp12+sp & yp3=yp1
;xp32=xp3+dx & yp32=yp3+dy
;;----------
;; right corner (top)
;xp4=xp22+sp & yp4=yp2
;xp42=xp4+dx & yp42=yp4+dy 
;;------------------------
;; left corner top                ; right corner top
;xpos(*,0) = [xp2,yp2,xp22,yp22] & xpos(*,1) = [xp4,yp4,xp42,yp42]
;; left corner bottom            ; right corner bottom
;xpos(*,2) = [xp1,yp1,xp12,yp12] & xpos(*,3) = [xp3,yp3,xp32,yp32]
;=============================================================
;======================================================
xpos      = fltarr(4,5)

dx=0.30 ;dx=0.22
dy=0.30 ;dy=0.40
sp=0.05
;---
; first left corner (bottom ) ; figure alone
xp0=0.05 & yp0=0.05
xp02=xp0+dx & yp02=yp0+dy
; left corner (bottom)
xp1=xp0 & yp1=yp02+sp
;xp1=0.15 & yp1=0.1
xp12=xp1+dx & yp12=yp1+dy
;-------
; left  corner (top)
xp2=xp1     & yp2=yp12+sp
xp22=xp2+dx & yp22=yp2+dy

;----------------
;right corner (bottom)
xp3=xp12+sp & yp3=yp1
xp32=xp3+dx & yp32=yp3+dy
;----------
; right corner (top)
xp4=xp22+sp & yp4=yp2
xp42=xp4+dx & yp42=yp4+dy 
;------------------------
;
; left corner top                ; right corner top
xpos(*,0) = [xp2,yp2,xp22,yp22] & xpos(*,1) = [xp4,yp4,xp42,yp42]
; left corner bottom            ; right corner bottom
xpos(*,2) = [xp1,yp1,xp12,yp12] & xpos(*,3) = [xp3,yp3,xp32,yp32]
; left corner bottom bottom
xpos(*,4)= [xp0,yp0,xp02,yp02]

;===================== x x x x x
;xpar      = ['!5qu','!5qv', '!5prate','!5T2m']
;xpar      = ['!5u10m','!5v10m', '!5prate','!5q2m']

npar      = n_elements(xpar)
;!p.multi = [0,2,2] & erase 
!p.multi = [0,2,3] & erase 
          
md=0 & xtit='!5ASIA - LIAM-1'; PC1
maplevels1 = indgen(10) * 0.1 + 0.1  
maplevels2 = indgen(10) * 0.1 - 1
    xmin1  = -1  & xmax1  = -0.4  & xmin2  = 0.4  & xmax2  = 1.
    low    = -1.0 & hig    = 1    & int    = 0.1
    color1 = [240,255,240,255]    & conlab = [1,0,1,0,1,0,1,0,1,0,1,0,1,0] 

; character size:
cs=1.5
; charthick
ctt=2.0
.run
for k = 0, npar - 1 do begin
tmp=cor(*,*,md,k)
;if(md eq 1) then tmp=-1*(tmp)
 pp=xpos(*,k)
  tvim,tmp,/noframe,barwidth=0.5,range=[0,1,2], $
  lcharsize=1.2,pos=pp,stitle=' ', labels='',colors=[255,255]
 contour,tmp,/overplot,/fill,levels= [xmin1,xmax1,xmin2,xmax2],$
    c_color=color1
 contour,tmp,/overplot,/noerase, levels=maplevels1,$
    C_THICK = 0.5, c_linestyle=0,color=2,c_labels=conlab
 contour,tmp,/overplot,/noerase, levels=maplevels2,$
    C_THICK = 0.5, c_linestyle=1,color=2,c_labels=conlab
 xyouts,mlon/2,mlat,xpar(k),alignment=0.5,charsize=1.25,color=2,charthick=1.5
 map_set,0,lonm,/cont,/usa,/cyl,limit=[lat1,lon1,lat2,lon2],/noerase, $
   charsize=1.3,/noborder,mlinethick=1.5,color=0,pos=xpos(*,k)
 if(k eq 0) then  $
    xyouts,lon1,lat2+5,xtit+'   (' +string(perc(md),format='(F4.1)')+'%)', $
      charsize=1.0,color=2,charthick=3
 axis,yaxis=0,ytype=0,ystyle=1,ycharsize=0.000001,yticklen=0.0000001
 axis,yaxis=1,ytype=0,ystyle=1,ycharsize=0.000001,yticklen=0.0000001
 axis,xaxis=0,xtype=0,xstyle=1,xcharsize=0.000001,xticklen=0.0000001
 axis,xaxis=1,xtype=0,xstyle=1,xcharsize=0.000001,xticklen=0.0000001
  axis,yaxis=0,ytype=0,ystyle=1,yrange=[lat1,lat2],charsize=cs,charthick=ctt
  axis,xaxis=0,xtype=0,xstyle=1,xrange=[lon1,lon2],charsize=cs,charthick=ctt
endfor 

;=====================================================
; PC2
erase

;!p.multi = [0,2,2] & erase 
!p.multi = [0,2,3] & erase 
          
md=1 & xtit='!5INDIA - LIAM2'; PC2
maplevels1 = indgen(10) * 0.1 + 0.1  
maplevels2 = indgen(10) * 0.1 - 1
    xmin1  = -1  & xmax1  = -0.3  & xmin2  = 0.3  & xmax2  = 1.
    low    = -1.0 & hig    = 1    & int    = 0.1
    color1 = [240,255,240,255]    & conlab = [1,0,1,0,1,0,1,0,1,0,1,0,1,0] 

for k = 0, npar - 1 do begin
tmp=cor(*,*,md,k)
;if(md eq 1) then tmp=-1*(tmp)
 pp=xpos(*,k)
  tvim,tmp,/noframe,barwidth=0.5,range=[0,1,2], $
  lcharsize=1.2,pos=pp,stitle=' ', labels='',colors=[255,255]
 contour,tmp,/overplot,/fill,levels= [xmin1,xmax1,xmin2,xmax2],$
    c_color=color1
 contour,tmp,/overplot,/noerase, levels=maplevels1,$
    C_THICK = 0.5, c_linestyle=0,color=2,c_labels=conlab
 contour,tmp,/overplot,/noerase, levels=maplevels2,$
    C_THICK = 0.5, c_linestyle=1,color=2,c_labels=conlab
 xyouts,mlon/2,mlat,xpar(k),alignment=0.5,charsize=1.25,color=2,charthick=1.5
 map_set,0,lonm,/cont,/usa,/cyl,limit=[lat1,lon1,lat2,lon2],/noerase, $
   charsize=1.3,/noborder,mlinethick=1.5,color=0,pos=xpos(*,k)
 if(k eq 0) then  $
    xyouts,lon1,lat2+5,xtit+'   (' +string(perc(md),format='(F4.1)')+'%)', $
      charsize=1.0,color=2,charthick=3
 axis,yaxis=0,ytype=0,ystyle=1,ycharsize=0.000001,yticklen=0.0000001
 axis,yaxis=1,ytype=0,ystyle=1,ycharsize=0.000001,yticklen=0.0000001
 axis,xaxis=0,xtype=0,xstyle=1,xcharsize=0.000001,xticklen=0.0000001
 axis,xaxis=1,xtype=0,xstyle=1,xcharsize=0.000001,xticklen=0.0000001
  axis,yaxis=0,ytype=0,ystyle=1,yrange=[lat1,lat2],charsize=cs,charthick=ctt
  axis,xaxis=0,xtype=0,xstyle=1,xrange=[lon1,lon2],charsize=cs,charthick=ctt
endfor 

  toggle
end 
;======================================================
;====================================================
; calculating duration, amplitude, onset, demise as in the iteractive
; program

;=============================================================================
;        smooth lisam1 and identify wet seasons
;
 window,xsize=1200,ysize=400 & loadct,5
;
.run /home/cjones/mjo/idl.for.subs/sub.persis.pro
;.run /home/cjones/mjo/idl.for.subs/filter.121.pro
.run /home/leila/idl_programs/subroutines/filter.121.annual.cycle.pro
;
; FOR PC1:
.run
  lisam1 = reform(pcs(*,0))
;  nti    = 15           ; <+++ +++ +++ 
 nti    = 3           ; <+++ +++ +++                    
; len     = 15
 len  = 5
  tsr    = lisam1

for m = 0, nti - 1 do begin
   tmp = smooth(tsr, len, /edge_truncate)
;   filter_121,tsr,mtot,nti,tmp
   tsr = tmp
endfor
        y  = intarr(mtot)
       idx = where(tsr gt 0)
    y(idx) = 1
   persistence,y, mtot, cluster,persis
mve,persis(where(persis ne 0))
end

!p.multi=[0,1,2] & erase
!p.multi(0) = 2
beg = 0  & ned = mtot-1

;----------------------------------
name=strcompress(string('lisam1_prec_tsr_pers_79.2013_len_',fix(len),'_nti_',fix(nti),'.'),/remove_all)
fim=strcompress(string(abs(fix(lat0)),'.N.',abs(fix(latf)),'N.lon',fix(lon0),'.',fix(lonf), '.eps'),/remove_all)

ffig=path2+name+fim

toggle,file=ffig,/landscape,/color

nada=tsr
nada(*)=0
;plot,cluster(beg:ned),yticklen=0.005,xrange=[0,ned-beg+1],xstyle=1,thick=1
plot,tsr,yticklen=0.005,xrange=[0,ned-beg+1],xstyle=1,thick=3,charsize=1.5,charthick=3,title='LIAM-1 smoothed'
oplot,nada,line=3,thick=2
!p.multi(0) = 1
plot,persis(beg:ned),yticklen=0.005,xrange=[0,ned-beg+1],xstyle=1,charsize=1.5,charthick=3,title='LIAM-1 persistence (days)',thick=3

toggle
;==========================================================================
;        crop to full wetseasons
;
mve,cluster(where(cluster ne 0))
;
   cook = cluster
 ; only for HS to remove the first and last cluster - wet season under-represented
;   idx = where(cluster eq 1)
;   cook(idx) = 0
;   idx = where(cluster eq max(cluster)) 
;   cook(idx) = 0
;****************************
   idx = where(cook ne 0)
    y  = intarr(mtot)
y(idx) = 1
persistence,y, mtot, wetseason,duration              ; retag 
;
!p.multi=[0,1,2] & erase
!p.multi(0) = 2
    beg = 0  & ned = mtot-1
plot,wetseason,yticklen=0.005,xrange=[0,mtot-1],xstyle=1
!p.multi(0) = 1
plot,duration,yticklen=0.005,xrange=[0,mtot-1],xstyle=1
mve,duration(where(duration ne 0))
;=============================================================================
;     monsoon amplitude: integral of positive values of PC1 
;     from onset to demise - NOTE: seasonal cycle is not removed
;
.run
   nseas       = max(wetseason)
   ampli1      = fltarr(mtot)               ; lisam1 amplitude
for m = 0, nseas - 1 do begin
   idx = where(wetseason eq m+1)
   ampli1(idx) = lisam1(idx)
endfor
end
;
!p.multi=[0,1,2] & erase 
!p.multi(0) = 2  & plot,ampli1,yticklen=0.005,xstyle=1

;==========================================================================
;     compute integrals of positive values of LISAM1
;     from onset to demise 
;
;
.run
   num        = max(wetseason)
   sams_ampli = fltarr(num)
for m = 0, num - 1 do begin
   idx  = where(wetseason eq m+1, pot)
   tsr  = ampli1(idx)
   x    = findgen(pot)
   idn  = where(tsr lt 0,cou)
   if(cou ne 0) then tsr(idn) = 0
   sams_ampli(m) = int_tabulated(x,tsr,/double)
endfor
end
;==========================================================================
;delvar,x,y,tsr,tmp,persis,beg,ned,cou,cook,cluster,idn,idx,m,num

;==========================================================================
;       print SAMS: onset,demise,duration,amplitude
;
rest=strcompress(string('len.',fix(len),'.nit.',fix(nti),'.dat'),/remove_all)
;filetxt=path2 + 'PC1.asia.amplitude.duration.prec.' + rest
;filetxt2=path2 + 'PC2.india.amplitude.duration.prec.'+ rest
fileamp1 = filetxt  + rest
fileamp2 = filetxt2 + rest

.run
  idx = where(cal_year eq 1999,ctot)               ;  create yearly calendar
  cli_day  = cal_day(idx)
  cli_mon  = cal_mon(idx)
  openw,1,fileamp1
  num = max(wetseason) & sep = strcompress(' |',/remove_all)
  sams_ons = intarr(3,num)                             ; day mon year
  sams_off = intarr(3,num)                             ;
  sams_dur = intarr(num)                               ;
for m = 0, num - 1 do begin
    idx           = where(wetseason eq m+1,dur)
    beg           = idx(0)
    ned           = idx(dur-1)
    id1 = where(cal_day(beg) eq cli_day and cal_mon(beg) eq cli_mon)
    id2 = where(cal_day(ned) eq cli_day and cal_mon(ned) eq cli_mon)
  print,m+1,cal_day(beg),cal_mon(beg),cal_year(beg),id1(0)+1,sep, $
          cal_day(ned),cal_mon(ned),cal_year(ned),id2(0)+1,dur, sep, $
          sams_ampli(m)
  sams_ons(0,m) = cal_day(beg) & sams_ons(1,m) =  cal_mon(beg)
  sams_ons(2,m) = cal_year(beg)
  sams_off(0,m) = cal_day(ned) & sams_off(1,m) =  cal_mon(ned)
  sams_off(2,m) = cal_year(ned)
  sams_dur(m)   = dur
  printf,1,format='(3I3,2I5,a3,2I3,3I5,a3,E15.5)', $
          m+1,cal_day(beg),cal_mon(beg),cal_year(beg),id1(0)+1,sep, $
          cal_day(ned),cal_mon(ned),cal_year(ned),id2(0)+1,dur, sep, $
          sams_ampli(m)
endfor
    close,1
end
;
;=============================================================================
;=============================================================================
; for PC2:


; FOR PC2:
.run

; lisam2 = reform(-1*pcs(*,1)) ; attention, using pc2 and multiplying by -1
  lisam2 = liam2
;  nti    = 15           ; <+++ +++ +++
;  nti    = 1                    
;  len     = 15

  tsr    = lisam2

for m = 0, nti - 1 do begin
   tmp = smooth(tsr, len, /edge_truncate)
;   filter_121,tsr,mtot,nti,tmp
   tsr = tmp
endfor
        y  = intarr(mtot)
       idx = where(tsr gt 0)
    y(idx) = 1
   persistence,y, mtot, cluster,persis
mve,persis(where(persis ne 0))
end

!p.multi=[0,1,2] & erase
!p.multi(0) = 2
beg = 0  & ned = mtot-1
nada=tsr
nada(*)=0
;plot,cluster(beg:ned),yticklen=0.005,xrange=[0,ned-beg+1],xstyle=1,thick=1

;----------------------------------

;----------------------------------
name=strcompress(string('liam2_79.2013.tsr_prec_pers_len_',fix(len),'_nti_',fix(nti),'.'),/remove_all)
fim=strcompress(string(abs(fix(lat0)),'.N.',abs(fix(latf)),'N.lon',fix(lon0),'.',fix(lonf), '.eps'),/remove_all)

ffig=path2+name+fim

toggle,file=ffig,/landscape,/color
nada=tsr
nada(*)=0
tit1=strcompress(string('LIAM-2 smoothed nit=',fix(nti),' len=',fix(len)))
tit2=strcompress(string('LIAM-2 persistence (days) nit',fix(nti),' len=',fix(len)))
;plot,cluster(beg:ned),yticklen=0.005,xrange=[0,ned-beg+1],xstyle=1,thick=1
plot,tsr,yticklen=0.005,xrange=[0,ned-beg+1],xstyle=1,thick=3,charsize=1.5,charthick=3,title=tit1
oplot,nada,line=3,thick=2
!p.multi(0) = 1
plot,persis(beg:ned),yticklen=0.005,xrange=[0,ned-beg+1],xstyle=1,charsize=1.5,charthick=3,title=tit2,thick=3

toggle
;========


;==========================================================================
;        crop to full wetseasons
;
mve,cluster(where(cluster ne 0))
;
   cook = cluster
 ; only for HS to remove the first and last cluster - wet season under-represented
;   idx = where(cluster eq 1)
;   cook(idx) = 0
;   idx = where(cluster eq max(cluster)) 
;   cook(idx) = 0
;****************************
   idx = where(cook ne 0)
    y  = intarr(mtot)
y(idx) = 1
persistence,y, mtot, wetseason2,duration2              ; retag 
;
!p.multi=[0,1,2] & erase
!p.multi(0) = 2
    beg = 0  & ned = mtot-1
plot,wetseason2,yticklen=0.005,xrange=[0,mtot-1],xstyle=1
!p.multi(0) = 1
plot,duration2,yticklen=0.005,xrange=[0,mtot-1],xstyle=1
mve,duration2(where(duration2 ne 0))
;=============================================================================
;     monsoon amplitude: integral of positive values of PC1 
;     from onset to demise - NOTE: seasonal cycle is not removed
;
.run
   nseas       = max(wetseason2)
   ampli2      = fltarr(mtot)               ; lisam1 amplitude
for m = 0, nseas - 1 do begin
   idx = where(wetseason2 eq m+1)
   ampli2(idx) = lisam2(idx)
endfor
end
;
!p.multi=[0,1,2] & erase 
!p.multi(0) = 2  & plot,ampli2,yticklen=0.005,xstyle=1

;==========================================================================
;     compute integrals of positive values 
;     from onset to demise 
;
;
.run
   num        = max(wetseason2)
   sams_ampli2 = fltarr(num)
for m = 0, num - 1 do begin
   idx  = where(wetseason2 eq m+1, pot)
   tsr  = ampli2(idx)
   x    = findgen(pot)
   idn  = where(tsr lt 0,cou)
   idp  = where(tsr gt 0,pcou)
   if(cou ne 0) then tsr(idn) = 0
;   print,m,cou,pcou
;    new version, adaptation for cases with persistence of one day
   if(pcou gt 1) then sams_ampli2(m) = int_tabulated(x,tsr,/double)
   if(pcou eq 0 and cou lt pot ) then sams_ampli2(m)=total(ampli2(idx))
   print,m,cou,pcou, sams_ampli2(m),pot
endfor
end

;==========================================================================
;delvar,x,y,tsr,tmp,persis,beg,ned,cou,cook,cluster,idn,idx,len,m,num

;==========================================================================
;       print : onset,demise,duration,amplitude for PC2
;
.run
  idx = where(cal_year eq 1999,ctot)               ;  create yearly calendar
  cli_day  = cal_day(idx)
  cli_mon  = cal_mon(idx)
  openw,1,fileamp2
  num = max(wetseason2) & sep = strcompress(' |',/remove_all)
  sams_ons2 = intarr(3,num)                             ; day mon year
  sams_off2 = intarr(3,num)                             ;
  sams_dur2 = intarr(num)                               ;
for m = 0, num - 1 do begin
    idx           = where(wetseason2 eq m+1,dur)
    beg           = idx(0)
    ned           = idx(dur-1)
    id1 = where(cal_day(beg) eq cli_day and cal_mon(beg) eq cli_mon)
    id2 = where(cal_day(ned) eq cli_day and cal_mon(ned) eq cli_mon)
  print,m+1,cal_day(beg),cal_mon(beg),cal_year(beg),id1(0)+1,sep, $
          cal_day(ned),cal_mon(ned),cal_year(ned),id2(0)+1,dur, sep, $
          sams_ampli2(m)
  sams_ons2(0,m) = cal_day(beg) & sams_ons2(1,m) =  cal_mon(beg)
  sams_ons2(2,m) = cal_year(beg)
  sams_off2(0,m) = cal_day(ned) & sams_off2(1,m) =  cal_mon(ned)
  sams_off2(2,m) = cal_year(ned)
  sams_dur2(m)   = dur
  printf,1,format='(3I3,2I5,a3,2I3,3I5,a3,E15.5)', $
          m+1,cal_day(beg),cal_mon(beg),cal_year(beg),id1(0)+1,sep, $
          cal_day(ned),cal_mon(ned),cal_year(ned),id2(0)+1,dur, sep, $
          sams_ampli2(m)
endfor
    close,1
end
;
;
;========== end PC2




;===== end procedure


;=============================================================================
;.............................................................................
; ;************************************************************
; saving data with indexes


openw,1,fffindx
.run
printf,1,'k  day1   mon1  year   liam1     liam2 '

for k=0,mtot-1 do begin

print,format='( 4(I5,1x), 2(F9.3,2x))',k,cal_day(k),cal_mon(k),cal_year(k),liam1(k),liam2(k)


printf,1,format='( 4(I5,1x), 2(F9.3,2x))',k,cal_day(k),cal_mon(k),cal_year(k),liam1(k),liam2(k)

endfor 
end

close,1


; save eofs; save eofs
;ff1=strcompress(string('lisam.trdEOF.prec.d.79.2013.CFSR.',abs(fix(lat0)),'N.',abs(fix(latf)),'N.lon',fix(lon0),'.',fix(lonf),'.IDL'),/remove_all)
; figure with eofs

;=============================================================================
delvar,beg,anom,avrg,clim,cluster,cook,data,idx,len,m,ned,nti,num,persis,sep
delvar,dur,qa2m,ta2m,tmp,tsmooth,tsr,u10m,v10m,y,yr1,totvar,modes,nfi,pcs,perc
delvar,res,i,j,l,aux,corre,spat,dim,prec,cor,x,pot,cou,idn

;==================================================
; endendendendendendendendend
;=============================================================================

restore,'/home/climate/wrf/hasia/trmm/TRMM3B42_hasia_1998_2010.IDL'
;restore,'/home/rkv-sbarc/arc/precipitation/trmm/TRMM3B42V7_dly00Z_1998.2013.HASIA.IDL'

tr_day=cal_day
tr_mon=cal_mon
tr_year=cal_year
tr_mtot=mtot
tr_mlat=mlat
tr_mlon=mlon
tr_yr1=yr1
tr_yr2=yr2
tr_rlat=rlat
tr_rlon=rlon
;
.run /home/leila/idl_programs/subroutines/Forest_Cannon/leap_year_map.pro
data=trmm
leap_year,data,tr_year,tr_mon,tr_day,tr_mtot,trmm


path2='/home/clivac/LISAM_results/Asia/india/CFSR_rain/2m10m/EOF_regular/1979_2013/'
; choose your variables here 
xpar      = ['!5Prec','!5u10m','!5v10m', '!5q2m','!5T2m']

fff=path2+ff1
restore,fff

id=where(cal_year ge tr_yr1 and cal_year le 2010,ntr) 
tr_lisam=liam1(id)
tr_lisam2=liam2(id)

idtr=where(tr_year ge tr_yr1 and tr_year le 2010,nnn)

ttrmm=trmm(*,*,idtr)

moyenvert,ttrmm,mtrmm,stdtrmm

atrmm=ttrmm


;====== calculating correlations with TRMM
tr_corr=fltarr(tr_mlon,tr_mlat,2)
tr_regr=fltarr(tr_mlon,tr_mlat,2)
.run
for i=0,tr_mlat-1 do begin
print,i
for j=0,tr_mlon-1 do begin
tsr=reform(ttrmm(j,i,*))
;   for k=0,tr_mtot-1 do begin
;    tsr(k)=tsr(k)-mtrmm(j,i)
; endfor 
res=correlate(tsr,tr_lisam)
tr_corr(j,i,0)=res(0)
res1 = REGRESS(tr_lisam, tsr, SIGMA=sigma, CONST=const, $
   MEASURE_ERRORS=measure_errors)
tr_regr(j,i,0)=res1(0)

res=correlate(tsr,tr_lisam2)
tr_corr(j,i,1)=res(0)

res1 = REGRESS(tr_lisam2, tsr, SIGMA=sigma, CONST=const, $
   MEASURE_ERRORS=measure_errors)
tr_regr(j,i,1)=res1(0)

endfor 
endfor 
end

!p.multi=[0,1,1]

ppp=[0.05,0.05,0.85,0.8]

;tr_lat1=min(tr_rlat)
;tr_lat2=max(tr_rlat)
;tr_lon1=min(tr_rlon)
;tr_lon2=max(tr_rlon)
;-------------------------
;tr_lat1=lat0
tr_lat1=-10
tr_lat2=latf
tr_lon1=lon0
;tr_lon2=lonf + 20.
tr_lon2=max(tr_rlon)

tmp=tr_rlat(0,*)

latid=where(tmp ge tr_lat1 and tmp le tr_lat2,nidla)

tmp=tr_rlon(*,0)
lonid=where(tmp ge tr_lon1 and tmp le tr_lon2,nidlo)
;------------------
; PLOTTING FANCE


; LOAD A COLOT TABLE
.run /home/fcannon/subroutines/color_tables.pro
ct,/anom
tvlct,r,g,b,/get
rgb = [[r],[g],[b]]

tmp1=tr_corr(lonid,latid,0)
rrlat=tr_rlat(0,latid)
rrlon=tr_rlon(lonid,0)

lati=rrlat(0)
lats=rrlat(-1)
loni=rrlon(0)
lons=rrlon(-1)
; LOAD DEM
restore,'/home/fcannon/elevation/highasia_DEM.IDL'

topo = congrid(hasiadem(3000:10200,0:4200,*),300,225,/interp,cubic=-.5)

; X&Y COORDINATES 
xxx = indgen(n_elements(tmp1(*,0)))*((lons-loni)/(n_elements(tmp1(*,0))-1))+loni
yyy = indgen(n_elements(tmp1(0,*)))*((lats-lati)/(n_elements(tmp1(0,*))-1))+lati

; OPEN WINDOW
w.close
; remember - this does not run with 'ei' command

w = window(window_title='LIAM Plots',Dimensions=[900,1100])
; SELECT MAXIMUM AND MINIMUM RANGE FOR YOUR DATA
ma=.7 & mi=-1*ma

; SET MAP
i1 = map('geographic',limit=[lati,loni,lats,lons+.1],/current,font_size=8,position=[.07,.67,.48,.9],/box_axes,grid_longitude=10,grid_latitude=10,linestyle=6,label_position=0)


i1['longitudes'].label_angle=0
; ADD Label
;Anomalies',font_size=10)
t1 = text(.085,.91,'a) Correlation TRMM and PC1',font_size=10)

mapp=tmp1
id=where(tmp1 ge -0.2 and tmp1 le 0.2,nid)
if(nid gt 0) then mapp(id)=0
; PLOT VARIABLE
c=image(mapp,xxx,yyy,/fill,grid_units='degrees',overplot=i1,rgb_table=rgb,min_value=mi,max_value=ma)
; OVERPLOT A FILLED CONTOUR OF PRECIPITABLE WATER ANOMALIES
;cccc=contour(tmp1,rrlon,rrlat,c_value=,grid_units='degrees',overplot=i1,c_value=indgen(20)*0.1-1.,c_color=indgen(1),c_label_show=1,c_thick=1)


; PLOT THE CONTINENTS
m1 = mapcontinents(color='black',/hires,thick=1)
; PLOT THE HIMALAYA (3000m CONTOUR)


;cc = contour(smooth(topo,1),indgen(300)*(60./(300-1))+45.,indgen(225)*(35./(225.-1)),rgb_table=0,overplot=i1,c_value=indgen(2)*3000+3000,/fill,c_color=indgen(1)+0,transparency=70,c_label_show=0,c_thick=2)
ccc = contour(smooth(topo,1),indgen(300)*(60./(300-1))+45.,indgen(225)*(35./(225.-1)),rgb_table=0,overplot=i1,c_value=indgen(2)*3000+3000,c_color=indgen(1)+0,c_label_show=0,c_thick=2)
; PLOT VARIABLE
;c=image(tmp1,xxx,yyy,/fill,grid_units='degrees',overplot=i1,rgb_table=rgb,min_value=mi,max_value=ma)


; ; TOP RIGHT MAP
tmp1=tr_corr(lonid,latid,1)
mapp=tmp1
id=where(tmp1 ge -0.2 and tmp1 le 0.2,nid)
if(nid gt 0) then mapp(id)=0

; SET MAP
i2 = map('geographic',limit=[lati,loni,lats,lons+.1],/current,font_size=8,position=[.5,.67,.91,.9],/box_axes,grid_longitude=10,grid_latitude=10,linestyle=6,label_position=0)
i2['longitudes'].label_angle=0
; ADD LABEL
t2 = text(.515,.91,'b) Correlation TRMM and PC-2 (LIMS)',font_size=10)
; PLOT VARIABLE
c=image(mapp,xxx,yyy,/fill,grid_units='degrees',overplot=i2,rgb_table=rgb,min_value=mi,max_value=ma)
;previous range = -8,8
; OVERPLOT THE WIND VECTORS
;vec = vector(mapu,mapv,rrlon,rrlat,overplot=i2,arrow_thick=2,head_proportional=1,length_scale=.75,head_size=.5,head_angle=0,x_subsample=2,y_subsample=2) ; the subsample keywords reduce the number of arrows for better visibility

; USE THE DATA TO CREATE A SAMPLE VECTOR FOR THE SCALE
;lll = legend(sample_magnitude=5,units=' m s!e-1!n',position=[rrlon(-1),rrlat(-1)],/data,vertical_alignment=0,target=vec,font_size=10,horizontal_alignment=0); the position might need to be changed according to your data's dimensions


; PLOT THE CONTINENTS
m1 = mapcontinents(color='black',/hires,thick=1)
; PLOT THE HIMALAYA (3000m CONTOUR)

; PLOT THE HIMALAYA (3000m CONTOUR)
;cc = contour(smooth(topo,1),indgen(300)*(60./(300-1))+45.,indgen(225)*(35./(225.-1)),rgb_table=0,overplot=i2,c_value=indgen(2)*3000+3000,/fill,c_color=indgen(1)+0,transparency=70,c_label_show=0,c_thick=2)
ccc = contour(smooth(topo,1),indgen(300)*(60./(300-1))+45.,indgen(225)*(35./(225.-1)),rgb_table=0,overplot=i2,c_value=indgen(2)*3000+3000,c_color=indgen(1)+0,c_label_show=0,c_thick=2)

; COLORBAR FOR TOP ROW
b.delete
;b = colorbar(position=[.905,.67,.92,.9],target=c,title='(kg m!e-2!n)',taper=1,border=1,font_size=10,orientation=1,textpos=1,tickinterval=2,tickdir=1,ticklen=2,minor=0)
b = colorbar(position=[.905,.67,.92,.9],target=c,title='Correlation',taper=1,border=1,font_size=10,orientation=1,textpos=1,tickinterval=0.1,tickdir=1,ticklen=2,minor=0)



;i4.save,'/home/clivac/LISAM_results/Asia/india/CFSR_rain/2m10m/EOF_regular/1979_2013/bogus_onset82_87.jpg',resolution=300
i2.save,'/home/clivac/LISAM_results/Asia/india/CFSR_rain/2m10m/EOF_regular/1979_2013/TRMM_LIMS_correlation_98_2010.jpg',resolution=300



;------------------------- old old old old old 
;;;;;;;;;;;;;;;;window,0

.run /home/leila/idl_programs/subroutines/loadb.pro
.run /home/leila/idl_programs/subroutines/loadidl.pro

loadct,5

toggle,file=ftrmm,/landscape,/color
erase
tmp1=tr_corr(lonid,latid,0)

tvim,tmp1,/scale,range=[-0.2,0.45],/noframe,barwidth=0.5, $
  lcharsize=1.2,pos=ppp,stitle=' ';, labels=''
;tlonm=(lon0+lonf)/2.
tlonm=(tr_lon1+tr_lon2)/2.
 map_set,0,tlonm,/cont,/cyl,limit=[tr_lat1,tr_lon1,tr_lat2,tr_lon2],/noerase, $
   charsize=1.3,/noborder,mlinethick=2,color=0,pos=ppp
 map_continents,/countries,/continents,/hires,color=0,pos=ppp

 axis,yaxis=0,ytype=0,ystyle=1,ycharsize=0.000001,yticklen=0.0000001
 axis,yaxis=1,ytype=0,ystyle=1,ycharsize=0.000001,yticklen=0.0000001
 axis,xaxis=0,xtype=0,xstyle=1,xcharsize=0.000001,xticklen=0.0000001
 axis,xaxis=1,xtype=0,xstyle=1,xcharsize=0.000001,xticklen=0.0000001
 axis,yaxis=0,ytype=0,ystyle=1,yrange=[tr_lat1,tr_lat2],charsize=1.5,charthick=4
 axis,xaxis=0,xtype=0,xstyle=1,xrange=[tr_lon1,tr_lon2],charsize=1.5,charthick=4

xyouts,tr_lon1,tr_lat2+2,'Correlation TRMM LIMS-1 (PC1)',font=5,charsize=2,charthick=2
 
erase
;;;; corr PC2 
;window,1
tmp2=tr_corr(lonid,latid,1)

tvim,tmp2,/scale,/noframe,range=[-0.2,0.45],barwidth=0.5, $
  lcharsize=1.2,pos=ppp,stitle=' ';, labels=''
tlonm=(tr_lon1+tr_lon2)/2.
 map_set,0,tlonm,/cont,/cyl,limit=[tr_lat1,tr_lon1,tr_lat2,tr_lon2],/noerase, $
   charsize=1.3,/noborder,mlinethick=2,color=0,pos=ppp
map_continents,/countries,/continents,/hires,color=0,pos=ppp

 axis,yaxis=0,ytype=0,ystyle=1,ycharsize=0.000001,yticklen=0.0000001
 axis,yaxis=1,ytype=0,ystyle=1,ycharsize=0.000001,yticklen=0.0000001
 axis,xaxis=0,xtype=0,xstyle=1,xcharsize=0.000001,xticklen=0.0000001
 axis,xaxis=1,xtype=0,xstyle=1,xcharsize=0.000001,xticklen=0.0000001
  axis,yaxis=0,ytype=0,ystyle=1,yrange=[tr_lat1,tr_lat2],charsize=1.5,charthick=4
  axis,xaxis=0,xtype=0,xstyle=1,xrange=[tr_lon1,tr_lon2],charsize=1.5,charthick=4

xyouts,tr_lon1,tr_lat2+2,'Correlation TRMM LIMS-2 (PC2)',font=5,charsize=2.,charthick=2


toggle
;-----------------

;------------------------------------------


exit
