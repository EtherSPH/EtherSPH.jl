'''
 # @ author: bcynuaa <bcynuaa@163.com> | callm1101 <Calm.Liu@outlook.com> | vox-1 <xur2006@163.com>
 # @ date: 2024/07/31 18:16:42
 # @ license: MIT
 # @ description:
 '''

import numpy as np
import matplotlib.pyplot as plt

plt.figure(facecolor="white", figsize=(8, 6))
plt.gca().set_aspect('equal', adjustable='box')
plt.axis('off')
#draw a pipe
pipe_width = 0.5
pipe_length = 2.0
plt.plot([0, 0], [0, pipe_width], 'k')
plt.plot([pipe_length, pipe_length], [0, pipe_width], 'k')
plt.plot([0, pipe_length], [0, 0], 'k')
plt.plot([0, pipe_length], [pipe_width, pipe_width], 'k')
plt.text(pipe_length/2, pipe_width/2, 'calculation domain', horizontalalignment='center', verticalalignment='center')
# draw buffer
buffer_length = 0.5
plt.plot([pipe_length, pipe_length + buffer_length], [0, 0], 'r')
plt.plot([pipe_length, pipe_length + buffer_length], [pipe_width, pipe_width], 'r')
plt.plot([pipe_length + buffer_length, pipe_length + buffer_length], [0, pipe_width], 'r')
plt.text(pipe_length + buffer_length/2, pipe_width/2, 'buffer', horizontalalignment='center', verticalalignment='center', color='r')
# draw an arrow from buffer to inlet
plt.hlines(pipe_width/2, pipe_length + buffer_length+0.1, pipe_length + buffer_length + 0.3, color='b')
plt.vlines(pipe_length + buffer_length + 0.3, pipe_width/2, pipe_width/2 + 0.4, color='b')
plt.hlines(pipe_width/2 + 0.4, -0.3, pipe_length + buffer_length + 0.3, color='b')
plt.vlines(-0.3, pipe_width/2, pipe_width/2 + 0.4, color='b')
plt.arrow(-0.3, pipe_width/2, 0.3, 0, color='b', head_width=0.05, head_length=0.05, length_includes_head=True)
plt.text(pipe_length/2, pipe_width/2 + 0.5, 'periodic', horizontalalignment='center', verticalalignment='center', color='b')
# inlet velocity shape
y = np.linspace(0, pipe_width, 100)
U = 0.8
plt.plot(U * y * (1 - y/pipe_width), y, color="g")
plt.plot(U * y * (1 - y/pipe_width) + buffer_length + pipe_length, y, color="g")
plt.text(0.15, pipe_width/2 - 0.1, r'$U_{\text{in}}$', color="g")
plt.text(0.15 + buffer_length + pipe_length, pipe_width/2 - 0.1, r'$U_{\text{in}}$', color="g")
# outlet velocity shape
y = np.linspace(0, pipe_width, 100)
U = 0.8
plt.plot(U * y * (1 - y/pipe_width) + np.random.uniform(-0.03, 0.03, len(y)) + pipe_length, y, color="g")
plt.text(pipe_length - 0.17, pipe_width/2 - 0.1, r'$U_{\text{out}}$', color="g")
# text
plt.text(pipe_length + buffer_length/2, -0.1, r'$U^*=(1-\xi)U_{\text{out}}+\xi U_{\text{in}}$', horizontalalignment='center', verticalalignment='center', color="r")
plt.savefig("example/cylinder/image/buffer_periodic_boundary_draft.png", bbox_inches='tight', dpi=300)