#!/bin/bash
#SBATCH --mem=4g
#SBATCH --time=7-00:00:00
#SBATCH --partition=long

sbatch run-1-200.sh

sleep 6h

sbatch run-201-400.sh

sleep 6h

sbatch run-401-600.sh

sleep 6h

sbatch run-601-800.sh

sleep 6h

sbatch run-801-1000.sh

sleep 6h

sbatch run-1001-1200.sh

sleep 6h

sbatch run-1201-1400.sh

sleep 6h

sbatch run-1401-1600.sh

sleep 6h

sbatch run-1601-1800.sh

sleep 6h

sbatch run-1801-2000.sh

sleep 6h

sbatch run-2001-2200.sh

sleep 6h

sbatch run-2201-2400.sh

sleep 6h

sbatch run-2401-2600.sh

sleep 6h

sbatch run-2601-2800.sh

sleep 6h

sbatch run-2801-3000.sh

sleep 6h

sbatch run-3001-3200.sh

sleep 6h

sbatch run-3201-3400.sh

sleep 6h

sbatch run-3401-3600.sh

sleep 6h

sbatch run-3601-3800.sh

sleep 6h

sbatch run-3801-4000.sh

sleep 6h

sbatch run-4001-4200.sh

sleep 6h

sbatch run-4201-4400.sh

sleep 6h

sbatch run-4401-4600.sh

sleep 6h

sbatch run-4601-4800.sh

sleep 6h

sbatch run-4801-5000.sh

sleep 6h

sbatch run-5001-5200.sh

sleep 6h

sbatch run-5201-5400.sh

sleep 6h

sbatch run-5401-5600.sh

sleep 6h

sbatch run-5601-5800.sh

sleep 6h

sbatch run-5801-6000.sh
