

## Usage 
**Developed in Unity 2019.3.11f1**

Open this directory in Unity or add the directory to your Unity Hub.

### GeneratorDemo.cs
Click on the Observatory object and enable only the GeneratorDemo.cs script. Press the play button then switch to the Scene view. The object will continue to be rendered at different angles and distances as dictated by the Update function. To stop the demo, simply click the play button again. 

*Note: You may need to adjust the perspective some to get a good view of the rendering process.*

### Generator.cs
Click on the Observatory object and enable only the Generator.cs script. Press the play button and wait as the dataset is generated. This process only takes about 4 minutes on an average laptop CPU. Once the script is complete you can navigate to the Dataset folder. The subdirectories will are named by the satellite images it contains. Each folder will contain 1000 224x224 pixel resolution PNG images of the satellite at different distances and orientations. To change the resolution or number of samples per class you must change the default arguments in the Start function.

*Note: Some satellites fail to correctly render all images, this problem is currently being debugged*
