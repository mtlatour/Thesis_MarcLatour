#!MC 1410

$!VARSET |epsilon| = 1e-6

$!VARSET |path| = "/Users/marclatour/Documents/Thesis/Thesis/Analysis/Tecplot_Output"

$!ACTIVEFIELDMAPS = [1]

$!VARSET |numPoints| = 100

$!VARSET |xMin| = |MINX|
$!VARSET |xMax| = |MAXX|

$!VARSET |yMin| = (|MINY|+|epsilon|)
$!VARSET |yMax| = (|MAXY|-|epsilon|)

$!VARSET |N| = 150

$!VARSET |dx| = ((|xMax| - |xMin|)/|N|)

$!LOOP |N|

        $!VARSET |xSlice| = ((|xMin| + (|LOOP|-1) * |dx|) + |epsilon|) 

        $!ACTIVEFIELDMAPS = [1]
        $!ExtendedCommand 
	       CommandProcessorID = 'Extract Precise Line'
	       Command = 'XSTART = |xSlice| YSTART = |yMin| ZSTART = 0 XEND = |xSlice| YEND = |yMax| ZEND = 0 NUMPTS = |numPoints| EXTRACTTHROUGHVOLUME = F EXTRACTTOFILE = F ' 

$!ENDLOOP

$!ACTIVEFIELDMAPS = [1-|NUMZONES|]

$!ExtendedCommand 
	CommandProcessorID = 'excsv'
	Command = 'VarNames:FrOp=1:ZnCount=|N|:ZnList=[2-|NUMZONES|]:AllVars:ValSep=",":FNAME="|path|/vertical-lines.csv"'
