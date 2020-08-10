import matplotlib.pyplot as plt
plt.close('all')
plt.subplots_adjust(0.1, 0.08, 0.96, 0.94, 0.2, 0.43)
plt.subplot(2, 2, 2 * ii + 1)
plt.imshow(20 * power_concrete,
extent=[times[0], times[-1], frequencies[0], frequencies[-1]],
aspect='auto', origin='lower', vmin=0., vmax=30.)
plt.xlabel('Time (s)')
plt.ylabel('Frequency (Hz)')
plt.title('Power (%s)' % title)
plt.colorbar()

plt.subplot(2, 2, 2 * ii + 2)
plt.imshow(phase_lock_concrete, extent=[times[0], times[-1], frequencies[0], frequencies[-1]], aspect='auto', origin='lower', vmin=0, vmax=0.7)
plt.xlabel('Time (s)')
plt.ylabel('Frequency (Hz)')
plt.title('Phase-lock (%s)' % title)
plt.colorbar()

plt.show()
