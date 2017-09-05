# sator
Analysis tool for the Arepo code.

Originally written by Thomas H. Greif and modified by Fernando Becerra.

### Running sator

Execute `./sator par.txt $function $base $SnapNumStart [$args]`

where `$function` can take the values:

0. Image
  * $args = $SnapNumEnd $ImgWidht
  
1. PhaseSpace
  * $args = $SnapNumEnd
  
2. Radial
  * $args = $SnapNumEnd $RadMin
  
3. M_BE
  * $args = $SnapNumEnd $MBEMin $ShiftX $ShiftY $ShiftZ
  
4. Halo
  
9. CutRegion
  * $args = $CutBase $CutSize
  
### Plotting

Execute `python python/$script $img_flag $base $SnapNumStart $SnapNumEnd` 

where `$script` can be:

0. image.py
  * $img_flag = 1
  
1. pspace.py
  * $img_flag = 0
  
2. radial.py
  * $img_flag = 0
  
3. mbe.py
  * $img_flag = 0
