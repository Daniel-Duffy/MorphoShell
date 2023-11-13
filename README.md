<p align="center">
  <a href="https://github.com/nschloe/dmsh"><img alt="dmsh" src="https://raw.githubusercontent.com/meshpro/dmsh/main/logo/logo-with-text.svg" width="50%"></a>
  <p align="center">The worst mesh generator you'll ever use.</p>
</p>

![plot](./banner_image.png)

Inspired by [distmesh](http://persson.berkeley.edu/distmesh/), dmsh can be slow,
requires a lot of memory, and isn't terribly robust either.


### Examples

#### Primitives

| <img alt="circle" src="https://raw.githubusercontent.com/meshpro/dmsh/main/plots/circle.svg" width="100%"> | <img alt="circle" src="https://raw.githubusercontent.com/meshpro/dmsh/main/plots/rectangle.svg" width="100%"> | <img alt="circle" src="https://raw.githubusercontent.com/meshpro/dmsh/main/plots/polygon.svg" width="100%"> |
| :--------------------------------------------------------------------------------------------------------: | :-----------------------------------------------------------------------------------------------------------: | :---------------------------------------------------------------------------------------------------------: |

```python
import dmsh
import meshio
import optimesh

geo = dmsh.Circle([0.0, 0.0], 1.0)
X, cells = dmsh.generate(geo, 0.1)

# optionally optimize the mesh
X, cells = optimesh.optimize_points_cells(X, cells, "CVT (full)", 1.0e-10, 100)

# visualize the mesh
dmsh.show(X, cells, geo)

# and write it to a file
meshio.Mesh(X, {"triangle": cells}).write("circle.vtk")
```

```python
import dmsh

geo = dmsh.Rectangle(-1.0, +2.0, -1.0, +1.0)
X, cells = dmsh.generate(geo, 0.1)
```

```python
import dmsh

geo = dmsh.Polygon(
    [
        [0.0, 0.0],
        [1.1, 0.0],
        [1.2, 0.5],
        [0.7, 0.6],
        [2.0, 1.0],
        [1.0, 2.0],
        [0.5, 1.5],
    ]
)
X, cells = dmsh.generate(geo, 0.1)
```

#### Combinations
