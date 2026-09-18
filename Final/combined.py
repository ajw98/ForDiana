#!/usr/bin/env python3

from PIL import Image
import os

# Desired final image size
final_width = 3900
final_height =1050  

# Spacing between images
spacing = 0

ncols = 4
nrows = 1



# Each tuple contains: (folder, filename)
image_files = [
    ("plt", "nfields_00000.png"),
    #("plt", "nfields_00005.png"),
    ("plt", "nfields_00078.png"),
    #("plt", "nfields_00085.png"),
    ("plt", "nfields_00090.png"),
    #("plt", "nfields_00095.png"),
    ("plt", "nfields_00100.png"),
    #("plt", "nfields_00103.png"),
    ("vel", "nvel_00000.png"),
    #("vel", "nvel_00005.png"),
    ("vel", "nvel_00078.png"),
    #("vel", "nvel_00085.png"),
     
    ("vel", "nvel_00090.png"),
    #("vel", "nvel_00095.png"),
    ("vel", "nvel_00100.png"),
    #("vel", "nvel_00103.png"),
]

# Size of each grid cell
cell_width = (final_width - spacing * (ncols - 1)) // ncols
cell_height = (final_height - spacing * (nrows - 1)) // nrows

# Create the final white image
combined = Image.new(
    "RGB",
    (final_width, final_height),
    "white"
)

for i, (folder, filename) in enumerate(image_files):

    image_path = os.path.join(folder, filename)

    image = Image.open(image_path).convert("RGB")

    # Resize the image to fit its grid cell
    image = image.resize(
        (cell_width, cell_height),
        Image.Resampling.LANCZOS
    )

    row = i // ncols
    col = i % ncols

    x = col * (cell_width + spacing)
    y = row * (cell_height + spacing)

    combined.paste(image, (x, y))

# Save the final image
combined.save(
    "combined_4x4.png",
    dpi=(600, 600)
)

print(f"Saved combined image: {final_width} × {final_height} pixels")
