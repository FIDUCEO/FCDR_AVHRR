#!/bin/sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.1.sh
echo python2.7 run_instrument.py NOAA06 >> run.1.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.2.sh
echo python2.7 run_instrument.py NOAA07 >> run.2.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.3.sh
echo python2.7 run_instrument.py NOAA08 >> run.3.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.4.sh
echo python2.7 run_instrument.py NOAA09 >> run.4.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.5.sh
echo python2.7 run_instrument.py NOAA10 >> run.5.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.6.sh
echo python2.7 run_instrument.py NOAA11 >> run.6.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.7.sh
echo python2.7 run_instrument.py NOAA12 >> run.7.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.8.sh
echo python2.7 run_instrument.py NOAA14 >> run.8.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.9.sh
echo python2.7 run_instrument.py NOAA15 >> run.9.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.10.sh
echo python2.7 run_instrument.py NOAA16 >> run.10.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.11.sh
echo python2.7 run_instrument.py NOAA17 >> run.11.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.12.sh
echo python2.7 run_instrument.py NOAA18 >> run.12.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.13.sh
echo python2.7 run_instrument.py NOAA19 >> run.13.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.14.sh
echo python2.7 run_instrument.py METOPA >> run.14.sh
bsub -q short-serial -W24:00 -oo run.1.log < run.1.sh
bsub -q short-serial -W24:00 -oo run.2.log < run.2.sh
bsub -q short-serial -W24:00 -oo run.3.log < run.3.sh
bsub -q short-serial -W24:00 -oo run.4.log < run.4.sh
bsub -q short-serial -W24:00 -oo run.5.log < run.5.sh
bsub -q short-serial -W24:00 -oo run.6.log < run.6.sh
bsub -q short-serial -W24:00 -oo run.7.log < run.7.sh
bsub -q short-serial -W24:00 -oo run.8.log < run.8.sh
bsub -q short-serial -W24:00 -oo run.9.log < run.9.sh
bsub -q short-serial -W24:00 -oo run.10.log < run.10.sh
bsub -q short-serial -W24:00 -oo run.11.log < run.11.sh
bsub -q short-serial -W24:00 -oo run.12.log < run.12.sh
bsub -q short-serial -W24:00 -oo run.13.log < run.13.sh
bsub -q short-serial -W24:00 -oo run.14.log < run.14.sh
