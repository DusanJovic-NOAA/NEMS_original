case `hostid` in

  0xac7c006)          MACHINE_ID=c ;; ### cirrus1
  0xac7c012)          MACHINE_ID=c ;; ### cirrus2

  0xac7a006)          MACHINE_ID=s ;; ### stratus1
  0xac7a012)          MACHINE_ID=s ;; ### stratus2

  0xffffffffbcc090b2) MACHINE_ID=g ;; ### gaea1
  0xffffffffbcc08fb2) MACHINE_ID=g ;; ### gaea2
  0xffffffffbcc08eb2) MACHINE_ID=g ;; ### gaea3
  0xffffffffbcc08db2) MACHINE_ID=g ;; ### gaea4

  ae0a2500)           MACHINE_ID=z ;; ### zeus1
  ae0a2600)           MACHINE_ID=z ;; ### zeus2
  ae0a2700)           MACHINE_ID=z ;; ### zeus3
  ae0a2800)           MACHINE_ID=z ;; ### zeus4
  ae0a2900)           MACHINE_ID=z ;; ### zeus5
  ae0a2a00)           MACHINE_ID=z ;; ### zeus6
  ae0a7d13)           MACHINE_ID=z ;; ### zeus7
  ae0a7e13)           MACHINE_ID=z ;; ### zeus8

  b20a057e)           MACHINE_ID=j ;; ### jet1
  b20a067e)           MACHINE_ID=j ;; ### jet2
  b20a077e)           MACHINE_ID=j ;; ### jet3
  b20a087e)           MACHINE_ID=j ;; ### jet4

esac

# --- for Zeus, find available account ID
if [ ${MACHINE_ID} = z ]; then
#  clear
#  ACCNR=null
#  for i in rm gm ada cmp gefs sepp omd
#  do
#    printf %s "Looking for accout ID: " $i " ..."
#    nr=`account_params 2>&1 | grep " "$i" " | wc -l`
#    if [ $nr -eq 1 ]; then
#      ACCNR=$i ; echo OK
#    else
#      echo
#    fi
#  done
#  if [ $ACCNR = null ]; then echo "Check your account ID"; exit ; fi
#  clear
ACCNR=gm
fi
