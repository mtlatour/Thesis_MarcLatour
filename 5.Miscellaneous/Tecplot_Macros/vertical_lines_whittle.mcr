#!MC 1410

$!VARSET |epsilon| = 1e-6

$!VARSET |path| = "/Users/marclatour/Documents/Thesis/Thesis/Analysis/Tecplot_Output"

$!ACTIVEFIELDMAPS = [1]

$!VARSET |numPoints| = 100

$!VARSET |xMin| = |MINX|
$!VARSET |xMax| = |MAXX|

$!VARSET |yMin| = (|MINY|+|epsilon|)
$!VARSET |yMax| = (|MAXY|-|epsilon|)

$!VARSET |N_1| = 140

$!VARSET |dx| = ((|xMax| - |xMin|)/|N_1|)

$!LOOP |N_1|

        $!VARSET |xSlice| = ((|xMin| + (|LOOP|-1) * |dx|) + |epsilon|) 

        $!ACTIVEFIELDMAPS = [1]
        $!ExtendedCommand 
	       CommandProcessorID = 'Extract Precise Line'
	       Command = 'XSTART = |xSlice| YSTART = |yMin| ZSTART = 0 XEND = |xSlice| YEND = |yMax| ZEND = 0 NUMPTS = |numPoints| EXTRACTTHROUGHVOLUME = F EXTRACTTOFILE = F ' 

$!ENDLOOP

$!ACTIVEFIELDMAPS = [2]

$!VARSET |numPoints| = 100

$!VARSET |xMin| = |MINX|
$!VARSET |xMax| = |MAXX|

$!VARSET |yMin| = (|MINY|+|epsilon|)
$!VARSET |yMax| = (|MAXY|-|epsilon|)

$!VARSET |N_2| = 60

$!VARSET |dx| = ((|xMax| - |xMin|)/|N_2|)

$!LOOP |N_2|

        $!VARSET |xSlice| = ((|xMin| + (|LOOP|-1) * |dx|) + |epsilon|) 

        $!ACTIVEFIELDMAPS = [2]
        $!ExtendedCommand 
	       CommandProcessorID = 'Extract Precise Line'
	       Command = 'XSTART = |xSlice| YSTART = |yMin| ZSTART = 0 XEND = |xSlice| YEND = |yMax| ZEND = 0 NUMPTS = |numPoints| EXTRACTTHROUGHVOLUME = F EXTRACTTOFILE = F ' 

$!ENDLOOP

$!ACTIVEFIELDMAPS = [2-|NUMZONES|]

$!VARSET |N| = (|N_1| + |N_2|)

$!ExtendedCommand 
	CommandProcessorID = 'excsv'
	Command = 'VarNames:FrOp=1:ZnCount=|N|:ZnList=[3-|NUMZONES|]:AllVars:ValSep=",":FNAME="|path|/vertical-lines.csv"'
