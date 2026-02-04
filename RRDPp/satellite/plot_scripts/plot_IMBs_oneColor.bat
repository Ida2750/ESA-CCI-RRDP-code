gmtset COLOR_FOREGROUND 183/9/25
gmtset COLOR_FOREGROUND 255/0/255
gmtset COLOR_BACKGROUND 255/0/255

set region=-90/60/90/60r
set proj=A-45/90/8c
set map_anot=g30/g10
set output=fig03

gmt psbasemap -R2001/2021/0/5 -JX6i/3i -Bpx2g1.0+l"Year" -By0.5g0.5+l"SIT [m]" -BWSne -K > %output%.ps
gawk -F"[&-]" "{print  $3+(($4+0.5)/12),$10}" IMB_info_table_Henriette.dat | gmt psxy -R -J -Sc0.15 -G89/89/89 -O -K >> %output%.ps
gawk -F"[- ]" "{print  $2+(($3+0.5)/12),$8}" SIMBA_initial_thickness.dat | gmt psxy -R -J -Sd0.20 -G89/89/89 -O >> %output%.ps
gmt psconvert %output%.ps -Tg -A