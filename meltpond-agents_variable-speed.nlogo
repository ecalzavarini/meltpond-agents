globals[
  ;mean-ice-height is defined via the GUI
  ;smooth-cycles   is defined via the GUI
  density-ratio
  ;melt-rate-ice
  ;melt-rate-pond
  water-height-transition
  ;seepage-rate

  gravity
  kinematic-viscosity
  horizontal-permeability
  albedo-sea-water
  albedo-ice

  ;; physical parameters for the model
  time-step
  space-step ; not used for the moment

  ;utility variables
  max1
  max2
]

breed[drops drop]

; internal variables on the patches:
; "ice" is ice height in [cm]
; "water" is the melt water in [cm]
; "melt" represents the water recently melted in [cm]
; "albedo" the fraction of radiation reflected by the surface
patches-own[ice water albedo]

drops-own[water-content]

;; we affect the global physical parameters that are relevant for the model
to startup
  set time-step  0.5  ; expressed in days
  set space-step 100.  ; the lateral size of a patch expressed in cm

  ;set melt-rate-ice 1.2  ;in cm of ice per day
  ;set melt-rate-pond 2.0 ;in cm of ice per day, nominal value, check the literature
  set water-height-transition 10 ;in cm
  set density-ratio 0.8 ;ratio between ice and water mass densities , for less porous ice can be 0.9

  ;ifelse seepage? [
  ;set seepage-rate 0.8 ; in cm / day
  ;][
  ;set seepage-rate 0.0 ; in cm / day
  ;]

  set gravity 9.81 * 100 * (86400 * 86400 ) ;in  cm / day^2
  set kinematic-viscosity 1.e-6 * (100 * 100) * 86400  ;in cm^2 / day
  set horizontal-permeability  3.e-9 * (100 * 100); in cm^2

  set albedo-sea-water 0.1
  set albedo-ice 0.9

  set max1 0
  set max2 0
end


; generate a random smooth topography
to setup-topography
  if clear-previous-plots? [clear-all]
  startup
   ask patches
  [
  set ice random-float (2 * mean-ice-height)
  set water 0
  set albedo compute-albedo

  set pcolor scale-color cyan ice 0 (2 * mean-ice-height)
  ]
  ifelse(smooth-with-radius?) [
    smooth-ice-with-radius][
    smooth-ice-with-cycles]
  reset-ticks
  tick ;this is to enable the plotting of the histogram of the initial configuration
  reset-ticks
end

; this is the function to smooth the random ice field
to smooth-ice-with-cycles
  repeat smooth-cycles [
    ask patches [
      let sum-ice-neighbors sum [ice] of neighbors4
      let mean-ice (sum-ice-neighbors + ice) / 5
      set ice mean-ice + random-float 0.1
      color-field
    ]
  ]
end

; this is the function to smooth the random ice field with a different procedure
to smooth-ice-with-radius
    ask patches [
      ;let sum-ice sum [ice] of patches in-radius smooth-radius
      ;let count-ice count patches in-radius smooth-radius
      ;let mean-ice (sum-ice  + ice) / ( count-ice + 1 )
      let mean-ice mean [ice] of patches in-radius smooth-radius
      set ice mean-ice
  ]
  ask patches[
      color-field
    ]
end



; coloring functions
to color-field
  ifelse water > 0 [
    set pcolor scale-color blue (2. * mean-ice-height * density-ratio - water) 0 (2. * mean-ice-height * density-ratio)
  ][
  if ice > 0 [set pcolor scale-color grey ice -10 (2 * mean-ice-height)]
  ;same result as the following 2 lines
  ;let ice-color 80 + ((89.9 - 80) / (2 * mean-ice-height)) * ice
  ;set pcolor ice-color
    if ice = 0 [set pcolor turquoise] ;blue - 3]
  ]
end



; melt ice
to melt-ice
  let actual-melted-volume 0
 if ice > 0 [
    ;; VERTICAL MELTING

    ;; 1) the following line implements conductive melting of ice
    if water = 0 [
      set actual-melted-volume  melt-rate-ice * time-step ;shall be corrected to include Stefan effect, see line below
      ;set actual-melted-volume  ( melt-rate-ice * time-step / ice * mean-ice-height )
    ]

    ;; 2) the following line implements the increased melt-rate for pondend ice : water enahnces melting (as in Luethje et al. paper)
    if water > 0 [
      ifelse water < water-height-transition[
      set actual-melted-volume (melt-rate-ice +  (melt-rate-pond - melt-rate-ice) * (water / water-height-transition)  )* time-step
      ][
      set actual-melted-volume melt-rate-pond * time-step
      ]
    ]

    ;; lateral melting
    if pond-lateral-melting?[set actual-melted-volume (actual-melted-volume + (melt-rate-pond * lateral-melting * time-step))]

    ;; Making melting happening
    if actual-melted-volume > ice [
      set actual-melted-volume ice  ;this is to avoid to melt more ice than what we have
      set water 0
    ]
    set ice (ice - actual-melted-volume ) ; melt has occurred
    sprout-drops 1 [
      set water-content (actual-melted-volume  * density-ratio)  ; the melted water is put into a pocket (a moving agent)
      ifelse pen-down? [pen-down][pen-up]
      hide-turtle
    ]
  ]
end

;; seepage of meltwater
to seepage
  let seepage-amount (seepage-rate * time-step) ; the seepage amount per patch

    ask drops-here [
      ifelse water-content > seepage-amount[
        set water-content water-content - seepage-amount
        set seepage-amount 0
      ][
        set seepage-amount (seepage-amount - water-content)
        set water-content 0
        die
      ]
    ]


  if seepage-amount > 0 and water > 0[
    ifelse water > seepage-amount[
      set water (water - seepage-amount)
    ][
      set water 0
    ]
  ]
end


;; move drops till conversion into water
;; flow assumes that the water produced by melt displaces instantaneously to position of minimum potential energy
to flow
  loop[
  let p min-one-of neighbors [ice + water]
   ifelse (ice + water) >  [ice + water] of p [
      move-to p
    ][
     set water (water + water-content) ; the pocket of melted water has reached a minimum height and it remains there
     if ice = 0 [set water 0]
     if not melt-ponds? [set water 0] ; remove all water if we are not interested in ponds
      die
  ]
  ]
end

; move drops till conversion into water
;; flow2 assumes that the spped is finite. It depends on the local gradient of "ice + water" and by the size of the patch. Needs improvement.
to flow2
  let p min-one-of neighbors [ice + water]
  ;let displacement (gravity / kinematic-viscosity) * (((ice + water) - [ice + water] of p) / space-step) * time-step
  ;let gradient (((ice + water) - [ice + water] of p) / space-step)
  let h atan space-step ((ice + water) - [ice + water] of p)
  let angle (90 - h) mod 360
  let horizontal-velocity  2 ;1000 * space-step  * sqrt (abs (gravity * (sin angle))  / (2 * space-step) )
  let displacement (horizontal-velocity * time-step)
   ifelse displacement >  distance p [
      move-to p
    ][
     set water (water + water-content) ; the pocket of melted water has reached a minimum height and it remains there
     if ice = 0 [set water 0]
     if not melt-ponds? [set water 0] ; remove all water if we are not interested in ponds
     die
  ]
end

to refreeze
  ask drops[
    set water (water + water-content) ; the pocket of melted water has reached a minimum height and it remains there
    if ice = 0 [set water 0]
    if not melt-ponds? [set water 0] ; remove all water if we are not interested in ponds
    die
  ]
  ask patches[
  set ice ice + (water / density-ratio)
  set water 0
  color-field
  ]
  tick
end

; a function for patches to account for the contribution to melting due to the presence of nearby ponds
to-report lateral-melting
  let ice-here ice
  let melt-volume 0
  ask neighbors4 with [ice < ice-here][
  ;ask neighbors with [ice < ice-here][
    ifelse (ice + water > ice-here)[
      set melt-volume (melt-volume + (ice-here - ice))
    ][
      set melt-volume (melt-volume + water)
    ]
  ]
  set melt-volume (melt-volume / space-step)
  report melt-volume
end



; this is the procedure for the loop over time
to melt-and-flow
  if pen-down? [cd] ; to clear previous drawings

  ask patches[
    ;; melt the ice
    melt-ice
    ;; seepage
    if seepage? [seepage]
  ]

  ask drops[
    ;; move water
    if water-flowing-mode = 1 [flow]
    if water-flowing-mode = 2 [flow2]
  ]

  ask patches[
    ;; recoloring the map
    color-field
    ;; estimate albedo
    set albedo compute-albedo
  ]

  if (count patches with [ice = 0] = count patches )[stop]
  tick
end



; reports used in plots
to-report compute-mean-ice
  report  mean [ice] of patches
end

to-report compute-mean-water
  report  mean [water] of patches
end

to-report compute-max-mean-water
  let water-depth mean [water] of patches
  if water-depth > max2 [ set max2 water-depth]
  report max2
end

to-report compute-std-ice
  let avg compute-mean-ice
  let var (mean [ice * ice] of patches)
  report sqrt (var - ( avg * avg ))
end

; we need to correct this, just a test
to-report compute-albedo
  let alpha albedo-sea-water  ; we start from total absorption (like if there is sea everywhere)
  if ice > 0 [
  ifelse water > 0 [
    let kappa (-1 / water-height-transition) * ln ( albedo-sea-water / albedo-ice )
    set alpha albedo-ice * exp (- water / water-height-transition) ; albedo coeff of pond
    ][
    set alpha albedo-ice ; albedo coeff of ice
    ]
  ]
  report alpha
end

to-report max-pond-coverage
  let value 100 * (count patches with [water > 0]) / count patches
  if value > max1 [set max1 value]
  report max1
end
@#$#@#$#@
GRAPHICS-WINDOW
230
48
742
561
-1
-1
4.0
1
10
1
1
1
0
1
1
1
0
125
0
125
1
1
1
ticks
30.0

BUTTON
20
211
216
244
NIL
setup-topography\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
20
67
215
100
mean-ice-height
mean-ice-height
0
400
93.0
1
1
cm
HORIZONTAL

SLIDER
20
175
215
208
smooth-cycles
smooth-cycles
0
40
3.0
1
1
NIL
HORIZONTAL

MONITOR
756
33
858
78
Average [cm]
compute-mean-ice
2
1
11

MONITOR
870
32
1011
77
Standard deviation [cm]
compute-std-ice
2
1
11

PLOT
1050
100
1326
308
Histogram of heights
[cm]
NIL
0.0
100.0
0.0
10.0
true
true
"set-plot-x-range 0 (mean-ice-height * 2)\n;set-plot-y-range 0 100\n;set-histogram-num-bars 7\n" ""
PENS
"ice" 1.0 0 -7500403 true "" "histogram [ice] of patches with [ice > 0]"
"water" 1.0 0 -13345367 true "" "histogram [water] of patches with [water > 0]"

TEXTBOX
757
10
907
28
Ice thickness
14
0.0
1

TEXTBOX
19
12
764
56
Arctic sea-ice evolution: melting, seepage and pond formation
18
0.0
1

TEXTBOX
19
45
215
75
1) create a random topography
12
105.0
1

BUTTON
546
583
655
616
NIL
melt-and-flow
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
756
100
1037
308
Average heights [cm]
time [day]
NIL
0.0
10.0
0.0
10.0
true
true
"set-plot-y-range 0 mean-ice-height" ""
PENS
"ice" 1.0 0 -7500403 true "" "plotxy (ticks * time-step) compute-mean-ice"
"water" 1.0 0 -13345367 true "" "plotxy (ticks * time-step) compute-mean-water"

TEXTBOX
1050
10
1200
28
Pond depth
14
0.0
1

MONITOR
1048
32
1143
77
Average [cm]
compute-mean-water
2
1
11

PLOT
756
316
1038
525
Relative coverage [%]
time [day]
NIL
0.0
10.0
0.0
100.0
true
true
"" ""
PENS
"water" 1.0 0 -13345367 true "" "plotxy (ticks * time-step) 100 * (count patches with [water > 0]) / count patches"
"ice" 1.0 0 -7500403 true "" "plotxy (ticks * time-step) 100 * (count patches with [ice > 0] ) / count patches"
"sea" 1.0 0 -14835848 true "" "plotxy (ticks * time-step) 100 * (count patches with [ice = 0] ) / count patches"

MONITOR
874
580
1024
625
Time to 100% melt [days]
ticks  * time-step
17
1
11

SLIDER
20
138
215
171
smooth-radius
smooth-radius
0
10
0.0
0.5
1
NIL
HORIZONTAL

SWITCH
20
102
215
135
smooth-with-radius?
smooth-with-radius?
1
1
-1000

TEXTBOX
235
565
385
593
3) run the simulation 
12
105.0
1

TEXTBOX
554
566
704
584
run step-by-step
11
15.0
1

TEXTBOX
669
565
893
593
run the complete simulation
11
15.0
1

BUTTON
662
582
845
615
NIL
melt-and-flow
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
405
583
533
616
pen-down?
pen-down?
1
1
-1000

SWITCH
21
330
217
363
melt-ponds?
melt-ponds?
0
1
-1000

SWITCH
22
448
219
481
seepage?
seepage?
0
1
-1000

PLOT
1056
316
1328
526
Albedo
time [day]
NIL
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"mean" 1.0 0 -2674135 true "" "plotxy (ticks * time-step) mean [albedo] of patches"
"ice" 1.0 0 -7500403 true "" "plotxy (ticks * time-step) albedo-ice"
"sea " 1.0 0 -14835848 true "" "plotxy (ticks * time-step) albedo-sea-water"

SWITCH
231
583
400
616
clear-previous-plots?
clear-previous-plots?
1
1
-1000

TEXTBOX
25
270
242
300
2) model features & parameters
12
105.0
1

MONITOR
1184
590
1327
627
diagnostoic count drops
count drops
17
1
9

CHOOSER
19
533
221
578
water-flowing-mode
water-flowing-mode
1 2
0

MONITOR
874
530
997
575
Max pond cover [%]
max-pond-coverage
2
1
11

SLIDER
22
405
218
438
melt-rate-pond
melt-rate-pond
0
10
2.0
0.1
1
cm/day
HORIZONTAL

SLIDER
21
291
217
324
melt-rate-ice
melt-rate-ice
0
10
1.2
0.1
1
cm /day
HORIZONTAL

SLIDER
21
486
220
519
seepage-rate
seepage-rate
0
10
0.8
0.1
1
H2O cm/day
HORIZONTAL

MONITOR
1154
32
1298
77
Max Average [cm]
compute-max-mean-water
2
1
11

SWITCH
22
367
218
400
pond-lateral-melting?
pond-lateral-melting?
1
1
-1000

@#$#@#$#@
# SUMMER MELTING OF ARCTIC SEA-ICE SHEETS 

## WHAT IS IT?

This is a highly idealized (toy?) model system for the evolution of melting dynamics of sea-ice sheets in the Arctic during the summer season.

### The model is based on the follwoing hypothesis:
1. Ice melt at a prescribed rate (described later on) and from this phase-change process melt water is generated. Note that water has a difference mass density as compared to ice
2. The water flows to the positions of local minima of potential energy, that is to say to positions where the sum of the ice thickness and the depth of melt water is minimum.
The water displacement process is assumed to occur over a time scale that is much smaller as compared to the time scale of the melting. This hypothesis allows to treat water down-slope movment as an istantaneous process.
3. When the ice thickness becomes zero, the melt water at the same location is taken out from the system. This feature of the model reproduces the discharge of fresh water at sea. 
4. seepage is accounted for

## HOW IT WORKS

### Ice-sheet surface generation
There are two functions that can be used to implement 

## HOW TO USE IT

1. Choose the characteristic of the ice sheet:
 * **mean-ice-height**:  _average thickness of the ice expressed in [cm]_ 
 * 


(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

The present is largely inspired from [1]


(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES


[1] M. Luethje, D. L. Feltham1, P. D. Taylor and M. G. Worster
_Modeling the summertime evolution of sea-ice melt ponds_
J. Geophys Res., **111**, C02001, 2006.

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
