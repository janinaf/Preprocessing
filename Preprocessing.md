# Documentation for Preprocessing.m

## Description
Organises raw behaviour data from SmoothWalk into individual laps and position bins, as well as align the timestamp with Intan.

## Requirements

SupportScripts:

loadintants.m

External:

NA 

## Usage
Preprocessing(mouseName,garrFile,intanFolder)

## Input
mouseName: name of mouse (string) | e.g. "MS14A"

garrFile: name of file where garr mat data is stored (string) | e.g. 'dis_vel_rew_2024-7-25_19_14_21_MS13B_REC_BASE.mat'

intanFolder: name of folder where binary intan data is stored (string) | e.g. 'W:\mouse_chemogenetics\MS14A\MS14A_240906_162727'

## Output

Output data automatically saved as "sess.mat" in the same directory as garrfile.

dataLap_all:
Table containing all behaviour data organised by lap and position bin
*Timestamps here NOT aligned to Intan

| Column | Variable |
| ------ | ------ |
| 1 | Lap Number |
| 2 | Position Bin |
| 3 | Timestamp |
| 4 | Track Position (Normalised) |
| 5 | Track Position (cm) |
| 6 | Instantaneous Velocity (cm/s) |
| 7 | Smoothed Velocity (cm/s) |
| 8 | Reward |
| 9 | Lick |
| 10 | Object Changes |
| 11 | Environment / Context |

dataLap:
Similar to dataLap_all, except
1) Only contains timestamps during movement
2) Timestamps aligned with Intan
3) Removed first and last incomplete lap
4) Removed fake laps
5) Lap numbers recomputed

dataBinned:
Time and velocity data organised into position bins
| Column | Variable |
| ------ | ------ |
| 1 | Lap Number |
| 2 | Position Bin |
| 3 | Start Time (s) |
| 4 | End Time (s) |
| 5 | Bin Duration (s) |
| 6 | Smoothed Velocity (cm/s) |

vel_binned:
Mean smoothed velocity per lap (row) and position bin (column)

nonzeroVel_idx:
Indices of dataLap_all where instantaneous velocity is not equal to zero

date:
Date extracted from garrfile name

offset:
Offset index from time alignment, use this to align LFP data

time_offset:
Time (s) of offset from time alignment, use this to align LFP data

reward:
| Column | Variable |
| ------ | ------ |
| 1 | Reward Timestamp (s) |
| 2 | Reward Position (Normalised) |

## Script workflow
1) Extract garr from MAT file generated from Smoothwalk, mainly need position information
2) Calculate and smooth velocity
3) Find indices of the start and end of each lap based on track position changes
4) Generate dataLap_all
5) Filter out timestamps with zero velocity (to be precise, filter using position change rather than instantaneous velocity)
6) Align timestamps of position information with Intan, and use Intan timestamps moving forward
7) Remove fake laps + first and last incomplete lap
8) Generate dataLap
9) Generate dataBinned
10) Extract reward signals from Intan (reward signals from Smoothwalk can be missing for some laps and more susceptible to imprecise timestamps)
11) Extract lick signals from Intan (*NOT DONE YET)
12) Save processed data as "sess.mat"
13) Save sanity check figures in "Preprocess_figures" subfolder (same directory as garrfile)