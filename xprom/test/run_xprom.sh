#!/bin/bash
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi
rm -f resultados* >& /dev/null
#awk 'BEGIN{printf "#Compound KinA KinB NetA NetB Interaction XC Clasic QA QB\n"}' > resultados.txt
echo "\begin{table}[!hbpt]" > resultados.txt
echo "  \centering" >> resultados.txt
echo "  \caption{IQA}" >> resultados.txt
echo "  \scalebox{0.8}{\begin{tabular}{l c c c c c c c c c}" >> resultados.txt
echo "    \hline \\\\" >> resultados.txt
echo "      \multicolumn{1}{l}{}" >> resultados.txt
echo "    & \multicolumn{1}{c}{\$K_{A}$}" >> resultados.txt
echo "    & \multicolumn{1}{c}{\$K_{B}$}" >> resultados.txt
echo "    & \multicolumn{1}{c}{\$E^{A}_{net}$}" >> resultados.txt
echo "    & \multicolumn{1}{c}{\$E^{B}_{net}$}" >> resultados.txt
echo "    & \multicolumn{1}{c}{\$E^{AB}_{Int}$}" >> resultados.txt
echo "    & \multicolumn{1}{c}{\$E^{AB}_{XC}$}" >> resultados.txt
echo "    & \multicolumn{1}{c}{\$E^{AB}_{Cl}$}" >> resultados.txt
echo "    & \multicolumn{1}{c}{\$Q_{A}$}" >> resultados.txt
echo "    & \multicolumn{1}{c}{\$Q_{B}$} \\\\" >> resultados.txt
echo "    \hline \hline \\\\" >> resultados.txt

for i in `ls -1 *.pmdout`; do
  name="${i%.pmdout}"
  $ECHO "  Running the parse for ${name} ... \c"
  xprom < options ${name}.pmdout resultados_${name}.txt >& /dev/null
  awk '/GROUPS 2 1 : Total Interaction Energy/ {inte=$9}
       /GROUPS 2 1 : XC Interaction Energy/ {xc=$9} 
       /GROUPS 2 1 : Classic Interaction Energy/ {cl=$9} 
       /GROUP  1   : Kinetic Energy/ {kina=$7}
       /GROUP  2   : Kinetic Energy/ {kinb=$7}
       /GROUP  1   : Net Energy/ {neta=$7}
       /GROUP  2   : Net Energy/ {netb=$7}
       /GROUP  1   : Total charge/ {qa=$7}
       /GROUP  2   : Total charge/ {qb=$7}
       END {printf "'${name}' & %10.6f & %10.6f & %10.6f & %10.6f & %10.6f & %10.6f & %10.6f & %10.6f & %10.6f \\\\ \n",kina,kinb,neta,netb,inte,xc,cl,qa,qb}' resultados_${name}.txt >> resultados.txt
  $ECHO "  Done"
done

echo "    \hline \\\\" >> resultados.txt
echo "  \end{tabular}}" >> resultados.txt
echo "  All in atomic units. cc-pVDZ basis set. All electrons correlated." >> resultados.txt
echo "  \label{iqa-table}" >> resultados.txt
echo "\end{table}" >> resultados.txt
