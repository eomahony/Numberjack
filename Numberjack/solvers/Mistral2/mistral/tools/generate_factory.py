#! /usr/bin/env python




map_var = [[['MaxWeight','FailureCountManager'],['wdeg','WDEG','Weighted Degree','WeightedDegree']],[['MinDomainOverWeight','FailureCountManager'],['dom/wdeg','dwdeg','DWDEG','Domain Over Weighted Degree','DomainOverWeightedDegree']],[['MaxWeight','ConflictCountManager'],['gwdeg','GWDEG','Global Weighted Degree','GlobalWeightedDegree']],[['MinDomainOverWeight','ConflictCountManager'],['dom/gwdeg','dgwdeg','DGWDEG','Domain Over Global Weighted Degree','DomainOverGlobalWeightedDegree']],[['MinWeight','ImpactManager'],['Impact','impact']],[['MinDomainTimesWeight','ImpactManager'],['IBS','ibs','Impact Based Search','impact based search']], [['MinDomainOverWeight','PruningCountManager'],['ABS','abs','activity based search','Activity Based Search','dom/activity']], [['MinDomain'], ['mindomain','MinDomain','min domain','minimum domain']], [['GenericNeighborDVO','SelfPlusAverage','MinDomainOverWeight', 'ConflictCountManager'], ['Neighbor','Neighbor','Neighbour','neighbour']], [['LexCombination< MinDomain, MaxDegree >'], ['mindom->maxdeg','MinDomainMaxDegree']], [['MinDomain'], ['mindomain','MinDomain','min domain','minimum domain','MinDomnMaxDeg']]]
           

#, [['LexCombination< MinDomain, MaxDegree >'], ['mindom->maxdeg','MinDomainMaxDegree']], [['MinDomain'], ['mindomain','MinDomain','min domain','minimum domain','MinDomnMaxDeg']]

           
map_val = [['AnyValue','No','no','Any','any'], ['MinValue','Lex','lex','minvalue','min value','MinValue','Min Value','lexicographic','Lexicographic'], ['MaxValue','AntiLex','antilex','maxvalue','max value','MaxValue','Max Value','antilexicographic','Antilexicographic'], ['HalfSplit','HalfSplit','halfsplit'],['RandomSplit','RandomSplit','RandSplit','randomsplit','randsplit'], ['RandomMinMax','RandomMinMax','randomminmax','RandMinMax','randminmax'], ['MinWeightValue', 'minweight', 'MinWeight'], ['Guided< MinValue >','MinVal+Guided','minval+guided'], ['Guided< MaxValue >','MaxVal+Guided','maxval+guided'], ['Guided< MinWeightValue >','MinWeight+Guided','minweight+guided'], ['Guided< MaxWeightValue >','MaxWeightVal+Guided','maxweightval+guided'], ['Guided< RandomMinMax >','Random+Guided','random+guided'], ['ConditionalOnSize< GuidedSplit< HalfSplit >, Guided< MinValue > >','Adpated','adapted']]


header = open('tools/ms.header', 'r')

footer = open('tools/ms.footer', 'r')

outfile = open('src/lib/mistral_solver.cpp', 'w')

for line in header:
    outfile.write(line)

for varo in map_var:
    outfile.write( '  if(' )
    VaroClass = varo[0]
    for vname in varo[1][:-1]:
        outfile.write( 'var_ordering == "'+vname+'" || ' )
    outfile.write( 'var_ordering == "'+varo[1][-1]+'" ) {\n' )
    for valo in map_val:
        ValoClass = valo[0]
        outfile.write( '    if( ')
        for vname in valo[1:-1]:
            outfile.write( 'branching == "'+vname+'" || ' )
        outfile.write( 'branching == "'+valo[-1]+'" ) {\n' )
        for r in range(1,6):
            if r == 1:
                outfile.write( '      if( randomness <= 1 ) {\n' )
            else:
                outfile.write( '      else if( randomness <= '+str(r)+' ) {\n' ) 
            if len(VaroClass) == 1:
                outfile.write( '        heu = new GenericHeuristic < GenericDVO < '+VaroClass[0]+' >, '+ValoClass+' > (this);\n' )    
            elif len(VaroClass) == 2:
                outfile.write( '        heu = new GenericHeuristic < GenericDVO < '+VaroClass[0]+', '+str(r)+', '+VaroClass[1]+' >, '+ValoClass+' > (this);\n' )          
            elif len(VaroClass) == 4:
                outfile.write( '        heu = new GenericHeuristic< '+VaroClass[0]+' < '+VaroClass[1]+', '+VaroClass[2]+', '+str(r)+', '+VaroClass[3]+' >, '+ValoClass+' > (this);\n' )
            outfile.write( '      }\n' )
        outfile.write( '    }\n' )
    outfile.write( '  }\n' )
        
            
for line in footer:
    outfile.write(line)

header.close()
footer.close()
outfile.close()
