##### Trying to Install MATLAB ######
# Cant get MATLAB to run

wget https://ssd.mathworks.com/supportfiles/downloads/R2019b/Release/4/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2019b_Update_4_glnxa64.zip
cd SSA*
mkdir MATLAB
mv ~/MAT* MAT*
cd MAT
unzip MAT*
rm *.zip
./install -mode silent -agreeToLicense yes
matlab 
    bash: matlab: command not found


##### Installing Octave 5.2.0 (Latest) #####

# Kind of installs 5.2.0 but running is weird: 
# https://flathub.org/apps/details/org.octave.Octave

sudo apt-get install flatpak
sudo apt install gnome-software-plugin-flatpak
sudo flatpak remote-add --if-not-exists flathub https://flathub.org/repo/flathub.flatpakrepo

# manually restart instance

sudo flatpak install flathub org.octave.Octave
flatpak run org.octave.Octave
    quit
alias octave="flatpak run org.octave.Octave"
octave
    plot(rand(5))
    F = getframe;
    error: getframe: not implemented for gnuplot graphics toolkit
    error: called from
        getframe at line 68 column 5

##### Installing Octave 4.0.3 ##### 

sudo apt-get install octave 

# Octave Kernel (allows you to dev octave in a notebook): 

pip install octave-kernel

# This will work but octave 4.0.3 doesn't support any of the camera functions


##### Force Octave 5.2.0 (apt install) #####

sudo apt-get install octave=5.2.0
Reading package lists... Done
Building dependency tree       
Reading state information... Done
E: Version '5.2.0' for 'octave' was not found