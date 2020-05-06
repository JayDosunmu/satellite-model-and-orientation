using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;
using System.Threading; // for hacky solution

// Author: Ryan James Walden

[RequireComponent(typeof(Camera))]
public class Generator : MonoBehaviour
{
    // Can probably remove these for OOP sake
    [SerializeField] private Transform _sunTransform;
    private Texture2D _screenshotTexture;
    private RenderTexture _renderTexture;
    private static string outputPath;
    private Camera _camera;

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
        int samples = 1000;

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

        _camera = GetComponent<Camera>();

         _sunTransform = GameObject.Find("Sun").transform;

        // Setup the observatory
        SetupTextures(resWidth, resHeight);

        StartCoroutine(GenerateSatellitesAndTakeImages(resWidth, resHeight, samples));
    }

    private IEnumerator GenerateSatellitesAndTakeImages(int resWidth, int resHeight, int samples)
    {

        // Get all GameObjects
        Object[] satelliteObjects = Resources.LoadAll("Satellites", typeof(GameObject));

        // Iterate through each object
        foreach (GameObject satellite in satelliteObjects)
        {
            // Update the console log
            Debug.Log("Rendering satellite: " + satellite.name);

            // Add a satellite to the scene
            GameObject mySatellite = Instantiate(satellite, Vector3.zero, Quaternion.identity);

            // Create directory for satellite if none exists
            outputPath = "Dataset/" + satellite.name;
            if (!Directory.Exists(outputPath)) Directory.CreateDirectory(outputPath);

            // Call the generation function
            //StartCoroutine(GenerateImages(mySatellite, resWidth, resHeight, samples));
            yield return StartCoroutine(GenerateImages(mySatellite, resWidth, resHeight, samples));

            // Remove the current Satellite from the scene, could make a corountine
            mySatellite.SetActive(false);
            Destroy(mySatellite);
            yield return null;
        }

        // Clean up camera
        CleanupObservatory();

        // Exit
        Debug.Log("Synthetic data generation complete");
        Application.Quit();
    }

    private void SetupTextures(int width, int height)
    {
        _renderTexture = new RenderTexture(width, height, 24);
        _screenshotTexture = new Texture2D(width, height, TextureFormat.RGB24, false);
        Debug.Log("Observatory setup complete");
    }

    Bounds CalculateBounds(GameObject obj)
    {
        Bounds b = new Bounds(obj.transform.position, Vector3.zero);
        Object[] rList = obj.GetComponentsInChildren(typeof(Renderer));
        foreach (Renderer r in rList)
        {
            b.Encapsulate(r.bounds);
        }
        return b;
    }

    void FocusObservatoryOnSatellite(GameObject satellite, float magnification = 1.0f)
    {
        Bounds b = CalculateBounds(satellite);
        Vector3 max = b.size;
        float radius = Mathf.Max(max.x, Mathf.Max(max.y, max.z));
        float dist = radius / (Mathf.Sin(_camera.fieldOfView * Mathf.Deg2Rad / 2f));
        _camera.transform.position = Random.onUnitSphere * dist * magnification + b.center;
        _camera.transform.LookAt(b.center);
    }


    private void CleanupObservatory()
    {
        _camera.targetTexture = null;
        RenderTexture.active = null;
        Destroy(_renderTexture);
        Debug.Log("Observatory cleanup complete");
    }


    //IEnumerator
    private IEnumerator GenerateImages(GameObject Satellite, int width, int height, int iterations)
    {
        // Grab the sun for transformations

        // Local for loop for now, want to distribute and parallelize eventually
        for (int i = 0; i < iterations; i++)
        {
            // We should only read the screen buffer after rendering is complete
            // Neccessary for IEnumerator 
            //yield return new WaitForEndOfFrame();

            // NOTE: All of the randoms need to have controlled variance for curriculum learning

            /* Randomly orient the sun
            GameObject.Find("Sun").transform.Rotate(new Vector3(Random.Range(-360.0f, 360.0f), 
                                                                Random.Range(-360.0f, 360.0f),
                                                                Random.Range(-360.0f, 360.0f)));
            */

            // Why do you rotate it here? lol
            // Randomly orient the satellite
            Satellite.transform.Rotate(new Vector3(Random.Range(0f, 360f),
                                                   Random.Range(0f, 360f),
                                                   Random.Range(0f, 360f)));

            // Focus the camera on the object
            FocusObservatoryOnSatellite(Satellite, Random.Range(0.75f, 2f));

            // Focus the sun on the object
            _sunTransform.rotation = _camera.transform.rotation;

            // Probably should change this
            // get main camera and manually render scene into rt
            _camera.targetTexture = _renderTexture;
            _camera.Render();

            // read pixels will read from the currently active render texture so make our offscreen 
            // render texture active and then read the pixels
            RenderTexture.active = _renderTexture;

            // Take a screenshot in 224x224 resolution
            // AsyncGPUReadback may be preffered to ReadPixels
            _screenshotTexture.ReadPixels(new Rect(0, 0, width, height), 0, 0);

            // Gets the correct size but unsure about location
            // Color[] img = ScreenShot.GetPixels(0, 0, resWidth, resHeight);

            // Save the screenshot to PNG in the corresponding satellite directory
            byte[] bytes = _screenshotTexture.EncodeToPNG();
            string filename = string.Format(outputPath + "/{0}x{1}_{2}.png", width, height, i);
            File.WriteAllBytes(filename, bytes);
            yield return null;
        }
    }
}
