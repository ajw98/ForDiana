#!/usr/bin/env python3

from PIL import Image

images = [
    Image.open("plt/nfields_00155.png"),
    Image.open("plt/nfields_00179.png"),
    Image.open("vel/nvel_00155.png"),
    Image.open("vel/nvel_00179.png"),
    ]


for i in range(2, 4):
    img = images[i]
    images[i] = img.crop((
        0,              # left
        200,            # top
        img.width,      # right
        img.height      # bottom
    ))

# Split into two rows
row1 = images[:2]
row2 = images[2:]

# Width of each column
col_widths = [
    max(row1[0].width, row2[0].width),
    max(row1[1].width, row2[1].width),
]

# Height of each row
row1_height = max(img.height for img in row1)
row2_height = max(img.height for img in row2)

total_width = sum(col_widths)
total_height = row1_height + row2_height

combined = Image.new(
    "RGB",
    (total_width, total_height),
    "white"
)

# First row
x = 0
for i, img in enumerate(row1):
    combined.paste(img, (x, 0))
    x += col_widths[i]

# Second row
x = 0
for i, img in enumerate(row2):
    combined.paste(img, (x, row1_height))
    x += col_widths[i]

combined.save("combined.png", dpi=(600, 600))
combined.save("combined.pdf", dpi=(600, 600))
