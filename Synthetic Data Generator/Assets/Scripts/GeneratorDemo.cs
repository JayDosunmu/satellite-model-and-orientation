using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;
using System.Threading; // for hacky solution

public class GeneratorDemo : MonoBehaviour
{

    private Texture2D ScreenShot;
    private RenderTexture rt;
    private GameObject mySatellite;
    private Transform mySun;

    // Start is called before the first frame update
    void Start()
    {
        Object[] satellites = Resources.LoadAll("Satellites", typeof(GameObject));
        mySatellite = Instantiate(satellites[0], Vector3.zero, Quaternion.identity) as GameObject;
        Thread.Sleep(250);
        mySun = GameObject.Find("Sun").transform;
        SetupObservatory(224,224);
    }

    // Update is called once per frame
    void Update()
    {
        Camera.main.targetTexture = rt;
        Camera.main.Render();
        RenderTexture.active = rt;
        mySatellite.transform.Rotate(new Vector3(Random.Range(0f, 360f), 
                                                   Random.Range(0f, 360f),
                                                   Random.Range(0f, 360f)));

        FocusObservatoryOnSatellite(Camera.main, mySatellite, Random.Range(0.4f, 1.5f));
        mySun.rotation = Camera.main.transform.rotation;
        Thread.Sleep(1000);
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
}
