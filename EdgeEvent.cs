using Unity.Mathematics;

namespace Briganti.StraightSkeletonGeneration
{
	public struct EdgeEvent
	{
		public int queueId;

		public EventType eventType;
		public float eventTime;
		public float2 eventPos;


		public void Reset()
		{
			// we never go back to being in the wavefront after being deleted from it
			if (eventType != EventType.NotInWavefront) eventType = EventType.None;
			eventTime = float.MaxValue;
		}

		public override string ToString()
		{
			if (eventType == EventType.NotInWavefront) return "";
			return $"{eventType} event at time {eventTime} and pos {eventPos}";
		}
	}
}