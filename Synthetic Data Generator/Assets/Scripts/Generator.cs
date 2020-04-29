using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;
using System.Threading; // for hacky solution

// Author: Ryan James Walden

public class Generator : MonoBehaviour
{
    // Can probably remove these for OOP sake
    private Texture2D ScreenShot;
    private RenderTexture rt;
    private static string outputPath;

    /* Fun things to consider:
        - Variation control
        - Caching
        - S3 / Cloud Storage
        - Lossless compression (crop to obj), save removed borders shape
        - Progress bar
        - More error handling and logging
        - Lighting variations
        - Artifiacts / celestial objects
        - Headless mode
    */


    void Start()
    {
        // Grab arguments
        string[] args = System.Environment.GetCommandLineArgs();

        // Setup default arguments
        int resWidth = 224;
        int resHeight = 224;
        int samples = 10;

        // Get args type
        var argsType = args.GetType();

        // Verify args type and args' element type
        if (argsType.IsArray)
        {
            // Verify all arguments are present
            if (args.Length != 3)
            {
                Debug.Log("Missing arguments");
            }
            else
            {
                // Setup arguments
                resWidth = System.Int32.Parse(args[0]);
                resHeight = System.Int32.Parse(args[1]);
                samples = System.Int32.Parse(args[2]);

                Debug.Log("Using custom arguments");
            }
        }
        else
        {
            Debug.Log("No arguments found");
        }
        
        // Setup the observatory
        SetupObservatory(resWidth, resHeight);

        // Get all GameObjects
        Object[] Satellites = Resources.LoadAll("Satellites", typeof(GameObject));

        // Iterate through each object
        foreach (GameObject Satellite in Satellites)
        {
            // Update the console log
            Debug.Log("Rendering satellite: " + Satellite.name);
            
            // Add a satellite to the scene
            GameObject mySatellite = Instantiate(Satellite, Vector3.zero, Quaternion.identity);

            // Create directory for satellite if none exists
            outputPath = "Dataset/" + Satellite.name;
            if (!Directory.Exists(outputPath)) Directory.CreateDirectory(outputPath);

            // Call the generation function
            //StartCoroutine(GenerateImages(mySatellite, resWidth, resHeight, samples));
            GenerateImages(mySatellite, resWidth, resHeight, samples);

            // Remove the current Satellite from the scene, could make a corountine
            mySatellite.SetActive(false);
            Destroy(mySatellite);
            //yield return new WaitForEndOfFrame();
        }
        
        // Clean up camera
        CleanupObservatory();
        
        // Exit
        Debug.Log("Synthetic data generation complete");
        Application.Quit();
    }


    private void SetupObservatory(int width, int height)
    {
        rt = new RenderTexture(width, height, 24);
        ScreenShot = new Texture2D(width, height, TextureFormat.RGB24, false);
        Debug.Log("Observatory setup complete");
    }


    Bounds CalculateBounds(GameObject obj) 
    {
        Bounds b = new Bounds(obj.transform.position, Vector3.zero);
        Object[] rList = obj.GetComponentsInChildren(typeof(Renderer));
        foreach (Renderer r in rList) {
            b.Encapsulate(r.bounds);
        }
        return b;
    }


    void FocusObservatoryOnSatellite(Camera c, GameObject satellite, float magnification = 1.0f) 
    {
        Bounds b = CalculateBounds(satellite);
        Vector3 max = b.size;
        float radius = Mathf.Max(max.x, Mathf.Max(max.y, max.z));
        float dist = radius / (Mathf.Sin(c.fieldOfView * Mathf.Deg2Rad / 2f));
        c.transform.position = Random.onUnitSphere * dist * magnification + b.center;
        c.transform.LookAt(b.center);
    }


    private void CleanupObservatory(){
        Camera.main.targetTexture = null;
        RenderTexture.active = null;
        Destroy(rt);
        Debug.Log("Observatory cleanup complete");
    }


    private static string ScreenShotName(int width, int height, int iter) 
    {
         return string.Format(outputPath + "/{0}x{1}_{2}.png", width, height, iter);
    }


    /* Apply Stuart's blur function to the images before saving them
    private void BlurImage(Array[])
    {

    }
    */

    //IEnumerator
    private void GenerateImages(GameObject Satellite, int width, int height, int itterations)
    {
        // Grab the sun for transformations
        Transform Sun = GameObject.Find("Sun").transform;
        
        // Local for loop for now, want to distribute and parallelize eventually
        for(int i=0;i<itterations;i++)
        { 
            // We should only read the screen buffer after rendering is complete
            // Neccessary for IEnumerator 
            //yield return new WaitForEndOfFrame();

            // Probably should change this
            // get main camera and manually render scene into rt
            Camera.main.targetTexture = rt;
            Camera.main.Render();

            // read pixels will read from the currently active render texture so make our offscreen 
            // render texture active and then read the pixels
            RenderTexture.active = rt;

            // NOTE: All of the randoms need to have controlled variance for curriculum learning

            /* Randomly orient the sun
            GameObject.Find("Sun").transform.Rotate(new Vector3(Random.Range(-360.0f, 360.0f), 
                                                                Random.Range(-360.0f, 360.0f),
                                                                Random.Range(-360.0f, 360.0f)));
            */

            // Randomly orient the satellite
            Satellite.transform.Rotate(new Vector3(Random.Range(0f, 360f), 
                                                   Random.Range(0f, 360f),
                                                   Random.Range(0f, 360f)));

            // Focus the camera on the object
            FocusObservatoryOnSatellite(Camera.main, Satellite, Random.Range(0.4f, 1.5f));

            // Focus the sun on the object
            Sun.rotation = Camera.main.transform.rotation;

            // Hacky solution
            Thread.Sleep(1000);

            // Take a screenshot in 224x224 resolution
            // AsyncGPUReadback may be preffered to ReadPixels
            ScreenShot.ReadPixels(new Rect(0, 0, width, height), 0, 0);
            
            // Gets the correct size but unsure about location
            // Color[] img = ScreenShot.GetPixels(0, 0, resWidth, resHeight);
            // data = BlurImage(img)

            // Save the screenshot to PNG in the corresponding satellite directory
            byte[] bytes = ScreenShot.EncodeToPNG();
            string filename = ScreenShotName(width, height, i);
            File.WriteAllBytes(filename, bytes);
            //yield return null;
        }
    }
}
