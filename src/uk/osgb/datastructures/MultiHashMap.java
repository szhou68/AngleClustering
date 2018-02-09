package uk.osgb.datastructures;
/**
 * this is a HashMap based simplified version of the multi-map data structure. 
 * 
 * At present it doesn't implement all methods in AbstractMap interface
 * 
 * 
 */

import java.util.AbstractSet;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

public class MultiHashMap<K, V> implements MultiMap<K, V>{
	
	private Map<K, HashSet<V>> map;
	//
	public MultiHashMap() {
		this(null);
	}
	// 
	public MultiHashMap(Map<K, V> copy) {
		map = new HashMap<K, HashSet<V>>();
		if (copy != null) {
			Iterator iter = copy.entrySet().iterator();
			while(iter.hasNext()) {
				Map.Entry<K, V> entry = (Map.Entry<K, V>)iter.next();
				put(entry.getKey(), entry.getValue());
			}
		}
	}
	
	public boolean containsKey(K key) {
		Collection values = (Collection)map.get(key);
		return ((values != null) && (values.size() != 0));
	}
	
	public boolean containsValue(V value) {
		Iterator iter = map.entrySet().iterator();
		boolean found = false;
		while (iter.hasNext()) {
			Map.Entry<K, HashSet<V>> entry = (Entry<K, HashSet<V>>) iter.next();
			Collection values = (Collection)entry.getValue();
			if (values.contains(value)) {
				found = true;
				break;
			}
		}
		return found;
	}
	//
	public Collection get(K key){
		return map.get(key);
	}
	// 
	public boolean put(K key, V value) {
		return getValues(key).add(value);
	}
	//
	public boolean putAll(K key, Collection values) {
		return getValues(key).addAll(values);
	}
	// remove all values associated with key
	public Collection remove(K key) {
		Collection original = getValues(key);
		map.remove(key);
		return original;
	}
	public boolean remove(K key, V value) {
		Collection values = (Collection)map.get(key);
		if (values == null) {
			return false;
		} else {
			return values.remove(value);
		}
	}
	//
	public String toString() {
		StringBuffer buff = new StringBuffer();
		buff.append("{");
		Iterator<K> keys = map.keySet().iterator();
		boolean first = true;
		while (keys.hasNext()) {
			if (first) {
				first = false;
			} else {
				buff.append(", ");
			}
			K key = keys.next();
			Collection values = getValues(key);
			buff.append("[" + key + ": " + values + "]");
		}
		buff.append("}");
		return buff.toString();
	}
	
	//
	public void clear() {
		map.clear();
	}
	//
	public int size(){
		return map.size();
	}
	//
	public Set<Map.Entry<K, V>> entrySet() {
		int size = 0;
		Iterator iterKeys = map.entrySet().iterator();
		while (iterKeys.hasNext()) {
			Map.Entry<K, HashSet<V>> entry = (Map.Entry<K, HashSet<V>>)iterKeys.next();
			Collection values = (Collection)entry.getValue();
			Iterator iterValues = values.iterator();
			while (iterValues.hasNext()) {
				size++;
				iterValues.next();
			}
		}
		final int finalSize = size;
		final Iterator entries = map.entrySet().iterator();
		return new AbstractSet() {
			int pos = 0;
			Map.Entry<K, HashSet<V>> entry;
			Iterator values;
			public Iterator iterator() {
				return new Iterator() {
					public void remove() {
						throw new UnsupportedOperationException();
					}
					public boolean hasNext() {
						return pos != finalSize;
					}
					public Map.Entry<K, V> next() {
						while(true) {
							if (entry == null) {
								entry = (Map.Entry<K, HashSet<V>>)entries.next();
								values = ((Collection)entry.getValue()).iterator();
							}
							K key = entry.getKey();
							if (values.hasNext()) {
								V value = (V) values.next();
								pos++;
								return new Entry(key, value);
							} else {
								entry = null;
							}
							
						}
					}
				};
			}
			public int size() {
				return finalSize;
			}
		};
	}	
	//
	private Collection getValues(K key) {
		HashSet<V> col = map.get(key);
		if (col == null) {
			col = new HashSet<V>();
			map.put(key, col);
		}
		return col;
	}
	//
	private static class Entry<K, V> implements Map.Entry<K, V> {
		K key;
		V value;
		Entry(K key, V value) {
			this.key = key;
			this.value = value;
		}
		public K getKey() {
			return key;
		}
		public V getValue() {
			return value;
		}
		public V setValue(V value) {
			V oldValue = this.value;
			this.value = value;
			return oldValue;
		}
		public boolean equals(Object o) {
			if (!(o instanceof Map.Entry)) {
				return false;
			} else {
				Map.Entry e = (Map.Entry)o;
				return (key==null ? e.getKey()==null : key.equals(e.getKey())) &&
				(value==null ? e.getValue()==null : value.equals(e.getValue()));
			}
		}
		public int hashCode() {
			return ((value==null) ? 0 : value.hashCode());
		}
		public String toString() {
			return key+"="+value;
		}
	}
	public void report(int maxEntries){
		if(map!=null){
			System.out.println("Number of entries: "+map.size());
			if(maxEntries <=0){
				maxEntries = map.size();
			}
			Iterator iterKeys = map.entrySet().iterator();
			int cnt = 0;
			while (iterKeys.hasNext()) {
				Map.Entry<K, HashSet<V>> entry = (Map.Entry<K, HashSet<V>>)iterKeys.next();
				Collection values = (Collection)entry.getValue();
				System.out.println(cnt+": "+values.size());
				cnt++;
				if(cnt >= maxEntries){
					return;
				}
			}
		}
	}
	@Override
	public Set<K> keySet() {
		// TODO Auto-generated method stub
		return map.keySet();
	}
}
