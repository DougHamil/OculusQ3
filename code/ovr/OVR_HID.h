#if !defined(_OVR_HID_H)
#define _OVR_HID_H

#include <linux/types.h>
#include <linux/input.h>
#include <linux/hidraw.h>

#ifndef HIDIOCSFEATURE
#error Update Kernel headers
#define HIDIOCSFEATURE(len)    _IOC(_IOC_WRITE|_IOC_READ, 'H', 0x06, len)
#define HIDIOCGFEATURE(len)    _IOC(_IOC_WRITE|_IOC_READ, 'H', 0x07, len)
#endif

#include <sys/ioctl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#include "OVR_Device.h"

#define OVR_VENDOR 0x2833
#define OVR_PRODUCT 0x0001

BOOLEAN getDeviceInfo( Device *dev );
BOOLEAN sendSensorScaleRange( Device *dev, const struct SensorScaleRange *r);
BOOLEAN sendSensorConfig(Device *dev, UInt8 flags, UInt8 packetInterval, UInt16 keepAliveIntervalMs);
BOOLEAN sendSensorKeepAlive(Device *dev);
BOOLEAN getSensorInfo( Device *dev );
BOOLEAN isRift( const char *path );
BOOLEAN openDevice(Device *dev, const char *path);
Device * openRiftHID( int nthDevice, Device *myDev );

#endif
