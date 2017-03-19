if [[ -n $1 ]]
then
  ext=$1
else
  ext=out
fi

echo -e "\nRsums"
grep "Rsum:" *.$ext

echo -e "\nSingle times"
grep "Threefry  Time single           round   2:" *.$ext

echo -e "\nMulti ss1 fix times"
grep "Threefry  Time multi ss1 fixlen round   2:" *.$ext

echo -e "\nMulti ss1 arg times"
grep "Threefry  Time multi ss1 arglen round   2:" *.$ext

echo -e "\nMulti all arg times"
grep "Threefry  Time multi all arglen round   2:" *.$ext

echo -e ""
