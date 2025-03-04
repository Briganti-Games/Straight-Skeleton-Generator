using Unity.Mathematics;

namespace Briganti.StraightSkeletonGeneration
{
	public enum SplitPoint
	{
		PrevVertex,
		Edge,
		NextVertex,
	}

	public struct EdgeEvent
	{
		public int queueId;

		public EventType eventType;
		public float eventTime;
		public float2 eventPos;

		public int reflexVertexIndex;

		public void Reset()
		{
			// we never go back to being in the wavefront after being deleted from it
			if (eventType != EventType.NotInWavefront) eventType = EventType.None;
			eventTime = float.MaxValue;
			reflexVertexIndex = -1;
		}

		public override string ToString()
		{
			if (eventType == EventType.NotInWavefront) return "";
			return $"{eventType} event at time {eventTime} and pos {eventPos}";
		}
	}
}