#!/bin/tcsh
set nstate=`grep 'Total R(velocity) tensor for State' $1 | wc -l`
echo $nstate > tensors.dat
echo $nstate > tensor-Rm.dat
echo $nstate > tensor-Rq.dat
grep ' nm ' $1 | cut -c39-47 >> tensors.dat
grep ' nm ' $1 | cut -c39-47 >> tensor-Rm.dat
grep ' nm ' $1 | cut -c39-47 >> tensor-Rq.dat
echo $nstate
if ( $nstate <= 9 ) then
  foreach x (`seq -w 1 $nstate`)
    grep 'Total R(velocity) tensor for State=          '$x -A4 $1 | tail -3 >> tensors.dat
    grep 'Rm(velocity) tensor for State=          '$x -A4 $1 | tail -3 >> tensor-Rm.dat
    grep 'Rq(velocity) tensor for State=          '$x -A4 $1 | tail -3 >> tensor-Rq.dat
  end
else if ( $nstate <= 99 ) then
  foreach x (`seq -w 1 9`)
    grep 'Total R(velocity) tensor for State=          '$x -A4 $1 | tail -3 >> tensors.dat
    grep 'Rm(velocity) tensor for State=          '$x -A4 $1 | tail -3 >> tensor-Rm.dat
    grep 'Rq(velocity) tensor for State=          '$x -A4 $1 | tail -3 >> tensor-Rq.dat
  end
  foreach x (`seq -w 10 $nstate`)
    grep 'Total R(velocity) tensor for State=         '$x -A4 $1 | tail -3 >> tensors.dat
    grep 'Rm(velocity) tensor for State=         '$x -A4 $1 | tail -3 >> tensor-Rm.dat
    grep 'Rq(velocity) tensor for State=         '$x -A4 $1 | tail -3 >> tensor-Rq.dat
  end

endif
#
