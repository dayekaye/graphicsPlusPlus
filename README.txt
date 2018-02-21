TO RUN:
navigate to directory of files
type "g++ -std=c++11 a7.cpp" (enter)
type "./a.out -f [filename.smf]" + any other options you want to include
open [filename.xps] to check solution

a6.cpp is the file that contains main and all my functions.  


Things I need to work on:
-3d clipping
-ironing out bugs in perspective view
-depth cueing

I started using c++11 for my assignment.  I wanted to get back in the swing of c++ since I know it is used a lot for graphics engineering. 
I am also running windows 10 pro with the ubuntu bash shell installed.

I used g++ to compile my assignment on my own computer.  this is what you should use as well.  I ran into no issues.

Had to do a lot to get my scanfilling working since I never had even 2d working properly.  
Implemented storing edges of a polygon so that I can pass that into the scan filling function.
I have back face culling working (which is apparent), and I have the depth buffer storing the correct z values. 
Unfortunately the depth buffer isn't shown since I couldn't finish the filling in time,
however if you look in the code you can see that I have this.


