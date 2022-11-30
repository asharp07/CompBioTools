### PyMOL Best Practices

These are some examples of some PyMOL visualizations I have made and example code I use to generate them. Using PyMOL to your full advantage is important with communicating your research. A low quality image can make it difficult for readers to follow your in-text explanations. 

<p align="center"><a href="https://imgur.com/IQCbgPu"><img src="https://i.imgur.com/IQCbgPu.png" width="400" /></a></p>

I always begin each PyMOL session with the following:

<p><code>remove hydro</code>
  
<code>remove solvent</code>
  
<code>remove resn CLA</code>
  
<code>remove resn POT</code>
  
<code>set bg_rgb, white</code>
  
<code>set ray_trace_mode, 1</code>
  
<code>set ray_trace_gain, 0.2</code>
  
<code>set antialias, 2 </code>
  
<code>set hash_max, 300 </code></p>

These ensure a crisp design and utilizing appropriate computational resources to generate a high quality image. I also like to make sure to use custom colors since some of the PyMOL defaults are less than exciting:

<code>set_color TEAL, [0, 102, 102]</code>

<code>color TEAL, selection</code>

This command colors the atoms in "selection" teal. For custom colors, I recommend using Adobe Color Wheel. I will be adding more example codes soon - including how to generate great looking membranes, protein-protein interfaces, and more! Below are a few examples of some of my favorite visualizations I've made during grad school. 

<a href="https://imgur.com/JJYocwP"><img src="https://i.imgur.com/JJYocwP.png" title="source: imgur.com" width="500" /></a>
<a href="https://imgur.com/vcbCOmK"><img src="https://i.imgur.com/vcbCOmK.png" title="source: imgur.com" width="500" /></a>
<a href="https://imgur.com/dGa0O45"><img src="https://i.imgur.com/dGa0O45.png" title="source: imgur.com" width="500" /></a>
<a href="https://imgur.com/CCUJUeC"><img src="https://i.imgur.com/CCUJUeC.png" title="source: imgur.com" width="500" /></a>
<a href="https://imgur.com/LkeBWGI"><img src="https://i.imgur.com/LkeBWGI.png" title="source: imgur.com" width="500" /></a>
<a href="https://imgur.com/h19kStS"><img src="https://i.imgur.com/h19kStS.png" title="source: imgur.com" width="500" /></a>
