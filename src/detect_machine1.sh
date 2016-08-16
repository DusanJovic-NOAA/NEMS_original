#!/bin/bash
#
# 06.11.16 E.Mirvis EMC/NCEP/NOAA - just to get a site
#
HOSTNAME='hostname -f'
export ACCNR=${ACCNR:-nems}

case `$HOSTNAME` in

  g10a1.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### gyre 1
  g10a2.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### gyre 2
  g14a1.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### gyre 3
  g14a2.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### gyre 4

  t10a1.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### tide 1
  t10a2.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### tide 2
  t14a1.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### tide 3
  t14a2.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### tide 4

  g20a1.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### gyre phase2
  g20a2.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### gyre phase2
  g20a3.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### gyre phase2
  g21a1.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### gyre phase2
  g21a2.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### gyre phase2
  g21a3.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### gyre phase2

  t20a1.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### tide phase2
  t20a2.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### tide phase2
  t20a3.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### tide phase2
  t21a1.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### tide phase2
  t21a2.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### tide phase2
  t21a3.ncep.noaa.gov)               MACHINE_ID=wcoss ;; ### tide phase2

  gaea1.ncrc.gov)                    MACHINE_ID=gaea ;; ### gaea1
  gaea2.ncrc.gov)                    MACHINE_ID=gaea ;; ### gaea2
  gaea3.ncrc.gov)                    MACHINE_ID=gaea ;; ### gaea3
  gaea4.ncrc.gov)                    MACHINE_ID=gaea ;; ### gaea4
  gaea5.ncrc.gov)                    MACHINE_ID=gaea ;; ### gaea5
  gaea6.ncrc.gov)                    MACHINE_ID=gaea ;; ### gaea6
  gaea7.ncrc.gov)                    MACHINE_ID=gaea ;; ### gaea7
  gaea8.ncrc.gov)                    MACHINE_ID=gaea ;; ### gaea8
  gaea9.ncrc.gov)                    MACHINE_ID=gaea ;; ### gaea9
  gaea10.ncrc.gov)                   MACHINE_ID=gaea ;; ### gaea10

  tfe01) MACHINE_ID=theia ;; ### theia01
  tfe02) MACHINE_ID=theia ;; ### theia02
  tfe03) MACHINE_ID=theia ;; ### theia03
  tfe04) MACHINE_ID=theia ;; ### theia04
  tfe05) MACHINE_ID=theia ;; ### theia05
  tfe06) MACHINE_ID=theia ;; ### theia06
  tfe07) MACHINE_ID=theia ;; ### theia07
  tfe08) MACHINE_ID=theia ;; ### theia08
  tfe09) MACHINE_ID=theia ;; ### theia09
  tfe10) MACHINE_ID=theia ;; ### theia10

esac

echo $MACHINE_ID    
