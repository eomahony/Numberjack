#echo "This script should be placed in the directory ~/entry_data"

cd /home/user/entry_data
echo "You should probably change the keyboard layout with \n dpkg-reconfigure keyboard-configuration (then restart)" 
echo "We assume that g++ and git are installed"

read -p "Press [Enter] key to start installation"

#apt-get install g++
#apt-get install git

git clone https://github.com/ehebrard/Mistral-2.0.git
cd Mistral-2.0
make clean 
make 
cd ..

cp Mistral-2.0/fz/mistral-fz fzn-exec

#echo "export PATH=\$DIR/entry_data/Mistral-2.0/fz:\$PATH" >> ../bin/challenge_env.sh
echo "PATH=/home/user/entry_data/Mistral-2.0/fz:$PATH" >> /home/user/.bashrc 

cp Mistral-2.0/fz/mznlib/* mzn-lib/
echo "Mistral-2.0 Installed.. We will performe lightweight tests."
read -p "Press [Enter] key to start test"
#nano fzn-exec

#
#OPTIONAL :  to install X
#apt-get install xfce4
#apt-get install gedit
#apt-get install Midori

#Test if mznlib is ok! 

cd

echo "exec-free black-hole 12 black-hole.mzn 12.dzn
exec-free black-hole 6 black-hole.mzn 6.dzn
exec-free on-call-rostering 10s-50d oc-roster.mzn 10s-50d.dzn" > entry_data/listfile


echo "Testing Global Constraints... [only cumulative]"

cd
wget homepages.laas.fr/msiala/minizinc2014/test--mzn-lib.sh 
chmod +x test--mzn-lib.sh 
./test--mzn-lib.sh 

ls 
echo "End"

