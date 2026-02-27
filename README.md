This is a C# implementation of the Straight Skeleton algorithm for Unity.

Straight Skeleton is an algorithm that takes a polygon (with holes) as an input and spits out a bunch of lines that resemble a roof of a building.

It basically turns this:

![A polygon with a hole in it](https://github.com/Briganti-Games/Straight-Skeleton-Generator/blob/master/Readme%20Images/polygon-input.jpg?raw=true)

Into this:

![A straight skeleton generated from a polygon](https://github.com/Briganti-Games/Straight-Skeleton-Generator/blob/master/Readme%20Images/straight-skeleton-output.jpg?raw=true)

This is very useful for procedurally generating roof structures.

There are several other implementations for straight skeletons out there, but all the ones I tried had several problems:
- Sometimes they were unstable for more complex polygons.
- Sometimes they were WAY too slow for realtime procedural generation.
- Sometimes they used too many resources or memory allocations.
- Sometimes they don't support polygons with holes.

This implementation does all of that. It has been battle-tested by thousands of users over the course of a year, because it's the core of the roof generation in [Dungeon Alchemist](https://www.dungeonalchemist.com/), the procedural map-generation app for roleplaying games that I'm working on. It's robust, works with very complex polygons with hundreds of vertices, is very fast and does not waste resources.

# Getting Started

Install this package using the Unity Package Manager via git URL.

First, go to **Window** > **Package Manager**.

Then, navigate to the following submenu:

![Screenshot of Unity Package Manager](https://github.com/Briganti-Games/Straight-Skeleton-Generator/blob/master/Readme%20Images/package-manager.png?raw=true)

Then, copy & paste the URL to this git package: https://github.com/Briganti-Games/Straight-Skeleton-Generator.git and press enter. This will start the import process.

Once imported, you can use the package.

This has been tested on Unity 2022.3, but should probably also work on earlier versions, as long as the Unity Mathematics package is supported, as this is the only dependency.

# How to use it?

# Implementation Details

