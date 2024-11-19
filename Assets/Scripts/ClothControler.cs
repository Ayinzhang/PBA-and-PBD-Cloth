using UnityEngine;
using UnityEngine.EventSystems;

public class ClothControler : MonoBehaviour, IDragHandler
{
    Transform clothTrans;

    void Start()
    {
        clothTrans = GameObject.Find("Cloth").transform;
    }

    public void OnDrag(PointerEventData eventData)
    {
        clothTrans.position += 0.001f * new Vector3(
            Mathf.Abs(clothTrans.position.x + 0.001f * eventData.delta.x) < 3? eventData.delta.x: 0, 
            Mathf.Abs(clothTrans.position.y + 0.001f * eventData.delta.y) < 1.5f? eventData.delta.y: 0, 
            0);
    }
}
