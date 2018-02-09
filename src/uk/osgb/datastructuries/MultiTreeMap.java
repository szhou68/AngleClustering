package uk.osgb.datastructures;
/**
 * this is a TreeMap-based simplified version of the multi-map data structure. 
 * 
 * At present, it doesn't implement all methods in AbstractMap interface 
 * 
 * 
 */

import java.util.AbstractSet;
import java.util.Collection;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

public class MultiTreeMap<K, V> implements MultiMap<K, V>{
	
	private TreeMap<K, TreeSet<V>> map;
	//
	public MultiTreeMap() {
		//this(null);
		map = new TreeMap<K, TreeSet<V>>();
	}
	public MultiTreeMap(Comparator<? super K> comparator){
		
	}
	//
	/*
	public MultiTreeMap(Map<K, V> copy) {
		map = new TreeMap<K, TreeSet<V>>();
		if (copy != null) {
			Iterator iter = copy.entrySet().iterator();
			while(iter.hasNext()) {
				Map.Entry<K, V> entry = (Map.Entry<K, V>)iter.next();
				put(entry.getKey(), entry.getValue());
			}
		}
	}
	*/
	public boolean containsKey(K key) {
		Collection values = (Collection)map.get(key);
		return ((values != null) && (values.size() != 0));
	}
	
	public boolean containsValue(V value) {
		Iterator iter = map.entrySet().iterator();
		boolean found = false;
		while (iter.hasNext()) {
			Map.Entry<K, TreeSet<V>> entry = (Entry<K, TreeSet<V>>) iter.next();
			Collection values = (Collection)entry.getValue();
			if (values.contains(value)) {
				found = true;
				break;
			}
		}
		return found;
	}
	public K ceilingKey(K key){
		return map.ceilingKey(key);
	}
	public K floorKey(K key){
		return map.floorKey(key);
	}
	public Collection<V> lastValue(){
		Map.Entry<K, TreeSet<V>> entry = map.lastEntry();
		return entry.getValue();
	}
	public Collection<V> firstValue(){
		Map.Entry<K, TreeSet<V>> entry = map.firstEntry();
		return entry.getValue();
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
	public void values(Collection<V> rlt){
		Set<Map.Entry<K, V>> entries = entrySet();
		for(Map.Entry<K, V> entry:entries){
			V val = entry.getValue();
			rlt.add(val);
		}
	}
	//
	public Set<Map.Entry<K, V>> entrySet() {
		int size = 0;
		Iterator iterKeys = map.entrySet().iterator();
		while (iterKeys.hasNext()) {
			Map.Entry<K, TreeSet<V>> entry = (Map.Entry<K, TreeSet<V>>)iterKeys.next();
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
			Map.Entry<K, TreeSet<V>> entry;
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
								entry = (Map.Entry<K, TreeSet<V>>)entries.next();
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
		TreeSet<V> col = map.get(key);
		if (col == null) {
			col = new TreeSet<V>();
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
				Map.Entry<K, TreeSet<V>> entry = (Map.Entry<K, TreeSet<V>>)iterKeys.next();
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
