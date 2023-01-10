more *.sdf | grep -A1 Smile | grep -v Smile > smile.dat
sed -i "s/--//g" smile.dat
grep -v "^$" smile.dat > smile_all.dat
