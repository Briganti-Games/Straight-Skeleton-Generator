using NUnit.Framework;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;

namespace Briganti.StraightSkeletons.Priority_Queue
{
	/// <summary>
	/// An implementation of a min-Priority Queue using a heap.  Has O(1) .Contains()!
	/// See https://github.com/BlueRaja/High-Speed-Priority-Queue-for-C-Sharp/wiki/Getting-Started for more information
	/// </summary>
	/// <typeparam name="T">The values in the queue.  Must extend the FastPriorityQueueNode class</typeparam>
	public class FastStructPriorityQueueTest
	{
		[Test]
		public static void Test()
		{
			FastStructPriorityQueue<int> queue = new FastStructPriorityQueue<int>(10);
			queue.Enqueue(1, 1, out int id1);
			queue.Enqueue(2, 2, out int id2);
			queue.Enqueue(3, 3, out int id3);

			Assert.IsTrue(queue.Dequeue() == 1);
			Assert.IsTrue(queue.Dequeue() == 2);
			Assert.IsTrue(queue.Dequeue() == 3);

			queue.Clear();
			queue.Enqueue(3, 3, out id3);
			queue.Enqueue(1, 1, out id1);
			queue.Enqueue(2, 2, out id2);

			Assert.IsTrue(queue.Dequeue() == 1);
			Assert.IsTrue(queue.Dequeue() == 2);
			Assert.IsTrue(queue.Dequeue() == 3);

			queue.Clear();
			queue.Enqueue(3, 3, out id3);
			queue.Enqueue(2, 2, out id2);
			queue.Enqueue(1, 1, out id1);

			Assert.IsTrue(queue.Dequeue() == 1);
			Assert.IsTrue(queue.Dequeue() == 2);
			Assert.IsTrue(queue.Dequeue() == 3);

			queue.Clear();
			queue.Enqueue(3, 3, out id3);
			queue.Enqueue(2, 2, out id2);
			queue.Enqueue(1, 1, out id1);
			queue.Remove(id2);

			Assert.IsTrue(queue.Dequeue() == 1);
			Assert.IsTrue(queue.Dequeue() == 3);
		}
	}
}